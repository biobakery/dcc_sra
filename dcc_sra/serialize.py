from os.path import basename
import xml.etree.ElementTree as ET
from xml.etree.ElementTree import SubElement as sub

from . import geo


def prep_subtype(p):
    return p._get_raw_doc()['meta']['subtype']

def eld(tagname, attrs={}, text=None, children=[]):
    ret = {"tagname": tagname}
    if attrs:
        ret['attrs'] = attrs
    if text:
        ret['text'] = text
    if children:
        ret['children'] = children
    return ret

def hier_sub(parent, tagname, attrs={}, text=None, children=[]):
    el = sub(parent, tagname, attrs)
    if text:
        el.text = text
    ret = [el]
    for kwargs in children:
        ret.append(hier_sub(el, **kwargs))
    return ret


def very_last(nested_list):
    last = nested_list[-1]
    while hasattr(last, "__iter__"):
        last = last[-1]
    return last


def spuid(obj):
    return eld("SPUID", attrs={"spuid_namespace": "hmp2"}, text=obj.id)


def flatten_list(l):
    def _flat(ls):
        for item in ls:
            if hasattr(item, "__iter__"):
                for subitem in _flat(item):
                    yield subitem
            else:
                yield item
    return list(_flat(l))


# Thanks, http://stackoverflow.com/a/4590052
def indent(elem, level=0):
    i = "\n" + level*"  "
    if len(elem):
        if not elem.text or not elem.text.strip():
            elem.text = i + "  "
        if not elem.tail or not elem.tail.strip():
            elem.tail = i
        for elem in elem:
            indent(elem, level+1)
        if not elem.tail or not elem.tail.strip():
            elem.tail = i
    else:
        if level and (not elem.tail or not elem.tail.strip()):
            elem.tail = i

def reg_text(t):
    return " ".join(t.decode('utf-8').split())

def reg_sample(s):
    s.mixs['lat_lon'] = " ".join(geo.cardinal(s.mixs['lat_lon']))
    return s


def _add_description(root, st):
    hier_sub(root, "Description", children=[
        eld("Comment", text="iHMP project "+st.name),
        eld("Organization", attrs={"role":"owner", "type":"institute"},
            children=[
                eld("Name", text="iHMP DCC"),
                eld("Contact", attrs={"email":"schwager@hsph.harvard.edu"},
                    children=[ eld("Name", children=[
                        eld("First", text="Randall"),
                        eld("Last",  text="Schwager")
                    ])
                ])
            ]
        )
    ])
    return root


def _add_bioproject(root, st):
    ret = hier_sub(root, "Action", children=[
        eld("AddData", attrs={"target_db":"BioProject"}, children=[
            eld("Data", attrs={"content_type":"xml"}, children=[
                eld("XmlContent", children=[
                    eld("Project", attrs={"schema_version":"2.0"})
                ])
            ]),
            eld("Identifier", children=[spuid(st)]),
        ])
    ])
    prj = flatten_list(ret)[-3]
    hier_sub(prj, "ProjectID", children=[spuid(st)])
    hier_sub(prj, "Descriptor", children=[
        eld("Title",    text="iHMP "+st.name),
        eld("Description", text=reg_text(st.description)),
        eld("Relevance", children=[ eld("Medical", text="Yes") ])
    ])
    pts_attrs = {"sample_scope":"eEnvironment"}
    hier_sub(prj, "ProjectType", children=[
        eld("ProjectTypeSubmission", attrs=pts_attrs, children=[
            eld("IntendedDataTypeSet", children=[
                eld("DataType", text="metagenome")
            ])
        ])
    ])
    return root


def _add_biosample(root, st, sample, prep):
    sample = reg_sample(sample)
    ret = hier_sub(root, "Action", children=[
        eld("AddData", attrs={"target_db":"BioSample"}, children=[
            eld("Data", attrs={"content_type":"xml"}, children=[
                eld("XmlContent", children=[
                    eld("BioSample", attrs={"schema_version":"2.0"})
                ])
            ]),
            eld("Identifier", children=[spuid(sample)])
        ])
    ])
    bs_node = flatten_list(ret)[-3]
    hier_sub(bs_node, "SampleId", children=[spuid(sample)])
    hier_sub(bs_node, "Descriptor", children=[
        eld("Title", text="%s %s sample"%(sample.mixs['env_package'],
                                          sample.mixs['body_product'])),
    ])
    hier_sub(bs_node, "Organism",
             attrs={"taxonomy_id": prep.ncbi_taxon_id},
             children=[eld("OrganismName", text="Metagenome")])
    hier_sub(bs_node, "Package", text="MIMS.me.human-associated.4.0")
    kv = lambda k, v: eld("Attribute", attrs={"attribute_name": k}, text=v)
    get = lambda v: sample.mixs.get(v, "missing")
    hier_sub(bs_node, "Attributes", children=[
        kv("env_biome", get("biome")),
        kv("collection_date", get("collection_date")),
        kv("env_feature", get("feature")),
        kv("env_material", get("material")),
        kv("geo_loc_name", get("geo_loc_name")),
        kv("host", "Homo sapiens"),
        kv("lat_lon", get("lat_lon"))
    ]+[kv(k, get(k)) for k in ("rel_to_oxygen", "samp_collect_device",
                               "samp_mat_process", "samp_size")
       if bool(sample.mixs.get(k, None))]
    )
    return root


def _add_sra(root, st, sample, prep, seq):
    kv = lambda k, v: eld("Attribute", attrs={"name": k}, text=v)
    strategy = "AMPLICON" if prep_subtype(prep) == "16s" else "WGS"
    mims_or_mimarks = prep.mimarks if prep_subtype(prep) == "16s" else prep.mims
    hier_sub(root, "Action", children=[
        eld("AddFiles", attrs={"target_db": "SRA"}, children=[
            eld("File", attrs={"file_path":basename(seq.urls[0])},
                children=[eld("DataType", text="generic-data")]),
            kv("instrument_model",seq.seq_model),
            kv("library_strategy",strategy),
            kv("library_source", "GENOMIC"),
            kv("library_selection", prep.lib_selection.upper()),
            kv("library_layout", "FRAGMENT"),
            kv("library_construction_protocol",
               reg_text(mims_or_mimarks['lib_const_meth'])),
            eld("AttributeRefId", attrs={"name": "BioProject"}, children=[
                eld("RefId", children=[spuid(st)])
            ]),
            eld("AttributeRefId", attrs={"name": "BioSample"}, children=[
                eld("RefId", children=[spuid(sample)])
            ]),
            eld("Identifier", children=[spuid(seq)])
        ])
    ])
    return root


def to_xml(st, samples):
    root = ET.Element('Submission')
    root = _add_description(root, st)
    root = _add_bioproject(root, st)
    for sample in samples:
        if not sample.prepseqs:
            continue
        for prep, seq in sample.prepseqs:
            root = _add_biosample(root, st, sample.sample, prep)
            root = _add_sra(root, st, sample.sample, prep, seq)

    return root

