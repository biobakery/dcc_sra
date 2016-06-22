import os
import time
import pipes
import socket
import string
import operator
from os.path import join, basename

import paramiko


last = operator.itemgetter(-1)

class SSHConnection(object):
    def __init__(self, user, host, keyfile, remote_path):
        self.remote_path = remote_path
        self.key = paramiko.RSAKey.from_private_key_file(keyfile)
        self.sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        self.sock.connect((host, 22))
        self.transport = paramiko.Transport(self.sock)
        self.transport.start_client()
        self.transport.auth_publickey(user, self.key)
        self.chan = self.transport.open_session()
        self.chan.get_pty()
        self.chan.invoke_shell()
        self._recvall()
        self._path_check()
        self.file_cache = self._build_file_cache()

    def _path_check(self):
        self.remote_path = self.remote_path.rstrip('/')
        head, tail = os.path.split(self.remote_path)
        ret = self.execute("ls -l "+head)
        fields = map(string.split, ret.strip().split('\r\n')[2:-1])
        names = dict(  zip(map(last, fields), fields) )
        if tail not in names:
            self.execute("mkdir --mode=0775 "+self.remote_path)
            return
        if tail in names and not names[tail][0].startswith("d"):
            self.execute("rm "+self.remote_path)
            self.execute("mkdir --mode=0775 "+self.remote_path)
            return


    def _build_file_cache(self):
        cache = dict()
        output = self.execute("ls -l "+self.remote_path).strip()
        for line in output.split('\r\n')[2:-1]:
            fields = line.split(None, 8)
            name = os.path.join(self.remote_path, fields[-1])
            try:
                size = int(fields[4])
            except ValueError:
                continue
            cache[name] = size
        return cache
            

    def _wait(self, tries=5):
        for i in range(tries):
            time.sleep(0.125)
            if self.chan.recv_ready():
                return True
        return False

        
    def _recvall(self):
        ret = str()
        while True:
            if self.chan.recv_ready():
                ret += self.chan.recv(1024)
            else:
                if not self._wait():
                    break
        return ret


    def execute(self, cmd, verbose=False):
        if not cmd.endswith("\n"):
            cmd += "\n"
        if verbose:
            print "sending `%s'"%(cmd)
        self.chan.send(cmd)
        return self._recvall()


    def fsize(self, fname):
        quoted = pipes.quote(join(self.remote_path, fname))
        output = self.execute("ls -l "+quoted).split("\n")[1].strip()
        try:
            return int(output.split(None, 8)[4])
        except ValueError:
            return None
        except IndexError:
            raise OSError(output)


    def uptodate(self, task, values):
        for fname in task.file_dep:
            if not os.path.exists(fname):
                return False
            remote_fname = join(self.remote_path, basename(fname))
            if remote_fname not in self.file_cache:
                return False
            return self.file_cache[remote_fname] == os.stat(fname).st_size

    def files(self):
        return self.execute("ls "+self.remote_path).strip().split('\r\n')[1:-1]
