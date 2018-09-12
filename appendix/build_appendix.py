import os
import re
from collections import OrderedDict
import markdown as md

PREAMBLE ='''\\newpage
\\section*{Appendix}
\\addcontentsline{toc}{chapter}{Appendix}
This section contains links to all the code and Jupyter notebooks used for generating figures. 
These can be found on GitHub at the base URL 
\\texttt{https://www.github.com/mparker2/mparker\_phd\_thesis/appendix}. \\\\
After navigating to this page, notebooks are organised by results chapter.
'''


def get_markdown_sections(fn):
    with open(fn) as f:
        text = f.readlines()
    m = md.Markdown()
    for prep in m.preprocessors.values():
        text = prep.run(text)
    tree = m.parser.parseDocument(text)
    items = iter(tree.getroot().getchildren())
    sections = OrderedDict()
    while True:
        try:
            h = next(items)
        except StopIteration:
            break
        if h.tag == 'h3':
            ol = next(items)
            assert ol.tag == 'ol'
            sections[h.text] = [c.getchildren()[0].text
                                for c in ol.getchildren()]
    return sections


def extract_section_notebook_urls(sections, base_url=None):
    sections_with_url = OrderedDict()
    for chapter_title, notebooks in sections.items():
        chapter_no, chapter_title = chapter_title.split(':')
        chapter_title, chapter_ref = chapter_title.split('[')
        chapter_ref = chapter_ref.replace(']', '')
        chapter_title = chapter_title.strip()
        subdir = chapter_no.lower().replace(' ', '_')
        if base_url:
            url = os.path.join(base_url, subdir)
        else:
            url = subdir
        notebooks_with_url = []
        for nb in notebooks:
            description, filename = nb.split(':')
            filename = filename.split(';')
            filename = [f.strip().replace('`', '') for f in filename]
            reference = filename[0].split('_', 1)[1]
            reference = re.sub(r'\.i?pyn?b?$', '', reference)
            filename = [(os.path.join(url, f), os.path.join(subdir, f)) for f in filename]
            filename = [(f[0].replace('_', '\\_'), f[1].replace('_', '\\_')) for f in filename]
            notebooks_with_url.append((description, filename, reference))
        sections_with_url[(chapter_title, chapter_ref)] = notebooks_with_url
    return sections_with_url

if __name__ == '__main__':
    appendix_sections = get_markdown_sections('appendix/README.md')
    appendix_sections = extract_section_notebook_urls(
        appendix_sections,
        base_url='https://www.github.com/mparker2/mparker_phd_thesis/appendix'
    )
    with open('appendix/appendix.tex', 'w') as o:
        o.write(PREAMBLE)
        o.write('\\begin{enumerate}\n')
        for i, (app_sect, notebooks) in enumerate(appendix_sections.items(), 2):
            chapter_title, chapter_ref = app_sect
            o.write('\\item{\n')
            o.write('\\textbf{{Chapter {}: {}}}\n'.format(i, chapter_title))
            o.write('\\begin{enumerate}[label*=\\arabic{*.}]\n')
            for description, filename, ref in notebooks:
                o.write('\\item{\n')
                o.write('{}:\\\\\n'.format(description))
                o.write('{} \\label{{{}}}\n'.format(
                    ';\\\\\n'.join(['\\texttt{{\\href{{{}}}{{{}}}}}'.format(*f) for f in filename]),
                    ref
                ))
                o.write('}\n')
            o.write('\\end{enumerate}\n')
            o.write('}\n')
        o.write('\\end{enumerate}\n')