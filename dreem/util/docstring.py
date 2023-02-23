
def docstring_to_dict_google(docstring):    
    docstring = docstring.replace( '    ', '\t')
    
    out = {}
    
    out['header'] = docstring.split('\n')[0]
    
    current_section = None
    for line in docstring.split('\n'):
        if line.replace('\t', '').replace(' ', '') == '':
            continue
        if line.startswith('\t\t\t') and current_section is not None:
            attr = line.strip()
            out[current_section][line.strip().split(' ')[0]] = ' '.join(line.strip().split(' ')[1:])
        elif line.startswith('\t\t'):
            current_section = line.strip()
            out[current_section] = {}
    
    return out


def dict_google_to_docstring(d):
    
    docstring = ''
    if 'header' in d.keys():
        docstring = d.pop('header')

    docstring += '\n'
    for key in d.keys():
        docstring += '\n\t' + key + '\n'
        for item in d[key]:
            docstring += '\t\t' + item + ' ' + d[key][item] + '\n'
    
    return docstring


def style_child_takes_over_parent(prnt_doc, child_doc): 

    docs = [prnt_doc, child_doc]
    for idx, doc in enumerate(docs):
        docs[idx] = docstring_to_dict_google(doc)
    prnt_doc, child_doc = docs
    
        
    for key in prnt_doc.keys():
        if key not in child_doc.keys():
            child_doc[key] = prnt_doc[key]
        else:
            if key == 'header':
                continue
            for item in prnt_doc[key]:
                if not item in child_doc[key].keys():
                    child_doc[key][item] = prnt_doc[key][item]

    return dict_google_to_docstring(child_doc)