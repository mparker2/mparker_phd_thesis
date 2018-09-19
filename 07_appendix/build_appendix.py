# hackiest script ever to turn README.md into appendix.tex

import os
import re
from collections import OrderedDict
import markdown as md

PREAMBLE ='''\\newpage
\\section*{Appendix}
\\addcontentsline{toc}{chapter}{Appendix}
This section contains links to all the code and Jupyter notebooks used for generating figures. 
These can be found on GitHub at the base URL 
\\texttt{https://www.github.com/mparker2/mparker\\_phd\\_thesis/tree/master/appendix}.\\\\
After navigating to this page, notebooks are organised by results chapter.
'''


def escape_underscores(url):
    if isinstance(url, str):
        return re.sub('_', '\\_', url)
    else:
        # assume its a list or tuple of strings
        return [escape_underscores(s) for s in url]


def clean_notebook_description(desc):
    return escape_underscores(desc.strip().strip(':'))


def get_label(notebook_name):
    _, filename = os.path.split(notebook_name)
    basename, _ = os.path.splitext(filename)
    basename = re.sub('^\d+\w?_', '', basename)
    return basename


def parse_markdown_to_xml_tree(fn):
    with open(fn) as f:
        text = f.readlines()
    m = md.Markdown()
    for prep in m.preprocessors.values():
        text = prep.run(text)
    root = m.parser.parseDocument(text).getroot()
    for treeprocessor in m.treeprocessors.values():
        new_root = treeprocessor.run(root)
        if new_root is not None:
            root = new_root
    return root


def get_markdown_sections(fn):
    root = parse_markdown_to_xml_tree(fn)
    items = iter(root.getchildren())
    sections = OrderedDict()
    while True:
        try:
            h = next(items)
        except StopIteration:
            break
        if h.tag == 'h3':
            section_name = escape_underscores(h.text)
            ol = next(items)
            assert ol.tag == 'ol'
            ol_items = []
            for list_item in ol.getchildren():
                list_item = list_item[0]
                assert list_item.tag == 'p'
                desc = clean_notebook_description(list_item.text)
                ahref = list_item.findall('a')
                notebooks = [a.text for a in ahref]
                urls = [a.get('href') for a in ahref]
                # only labels don't want underscores escaped
                label = get_label(notebooks[0])
                urls = escape_underscores(zip(urls, notebooks))
                ol_items.append((desc, urls, label))
            sections[section_name] = ol_items
    return sections


if __name__ == '__main__':
    appendix_dir = os.path.split(os.path.realpath(__file__))[0]
    fn = os.path.join(appendix_dir, 'README.md')
    appendix_sections = get_markdown_sections(fn)
    with open(os.path.join(appendix_dir, 'appendix.tex'), 'w') as o:
        o.write(PREAMBLE)
        o.write('\\begin{enumerate}\n')
        for app_sect, notebooks in appendix_sections.items():
            chapter_title = app_sect
            o.write('\\item{\n')
            o.write('\\textbf{{{}}}\n'.format(chapter_title))
            o.write('\\begin{enumerate}[label*=\\arabic*]\n')
            for description, filename, ref in notebooks:
                o.write('\\item{\n')
                o.write('{}:\\\\\n'.format(description))
                urls = ';\\\\\n'.join(
                    ['\\texttt{{\\href{{{}}}{{{}}}}}'.format(*f) for f in filename]
                )
                o.write('{} \\label{{{}}}\n'.format(urls, ref))
                o.write('}\n')
            o.write('\\end{enumerate}\n')
            o.write('}\n')
        o.write('\\end{enumerate}\n')