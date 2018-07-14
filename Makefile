clean :
	rm -f */tmp.md
	rm -f */figures/*.png
	rm -f title_pages/acronyms.sorted.csv

pdf :
	# replace figure urls to include chapter directory and png
	# also add bold text at end for pandoc-shortcaption
	for CHAPTER_TEXT in $$(ls */text.md); do \
		CHAPTER_NAME="$${CHAPTER_TEXT%%/text.md}"; \
		CHAPTER_TEMP="$${CHAPTER_NAME}/tmp.md"; \
		sed -e "s/figures\/\(.*\).svg/$${CHAPTER_NAME}\/figures\/\1.png/" $$CHAPTER_TEXT | \
		sed -e "/^!\[/s/%/\\\\\\\\%/g" -e "/^!\[/s/&/\\\\\\\\&/g" | \
		perl -pe 's|(!\[\*\*(.*?)\*\*.*png)|\1 "\2"|' > $$CHAPTER_TEMP; \
	done
	# create png versions of svgs
	for IMG in $$(ls */figures/*.svg); do \
		inkscape -d 300 -e $${IMG%%.svg}.png $$IMG; \
	done
	sort -k1,1 -t',' title_pages/acronyms.csv > title_pages/acronyms.sorted.csv
	# pandoc doesn't accept markdown include-before-body with tex template
	# SO we have to render the abstract first
	pandoc \
		-o "title_pages/tmp.tex" \
		--write=latex \
		--read=markdown \
		title_pages/acknowledgements.md \
		title_pages/abstract.md
	pandoc \
		-o "mparker_thesis_draft_`date +%d%m%y`.pdf" \
		--template "style/template.tex" \
		-V fontsize=12pt \
		-V papersize=a4paper \
		-V documentclass=report \
		--pdf-engine=xelatex \
		--filter ./pandoc-shortcaption/pandoc_shortcaption/__main__.py \
		--read=markdown \
		--include-before-body title_pages/title_pages.tex \
		--include-before-body title_pages/tmp.tex \
		--bibliography references/references.bib \
		--csl references/cell.csl \
		title_pages/contents.tex \
		introduction/tmp.md \
		chapter_3/tmp.md \
		chapter_4/tmp.md \
		chapter_5/tmp.md \
		chapter_6/tmp.md \
		discussion/tmp.md \
		references/references.tex ;
	rm */tmp.md
	rm */tmp.tex
	rm */figures/*.png
	rm title_pages/acronyms.sorted.csv