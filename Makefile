clean :
	rm */tmp.md
	rm */figures/*.png
full :
	# replace figure urls to include chapter directory and png
	for CHAPTER_TEXT in $$(ls */text.md); do \
		CHAPTER_NAME="$${CHAPTER_TEXT%%/text.md}"; \
		CHAPTER_TEMP="$${CHAPTER_NAME}/tmp.md"; \
		echo $$CHAPTER_TEMP; \
		sed -e "s/figures\/\(.*\).svg/$${CHAPTER_NAME}\/figures\/\1.png/" $$CHAPTER_TEXT > $$CHAPTER_TEMP; \
	done
	# create png versions of svgs
	for IMG in $$(ls */figures/*.svg); do \
		inkscape -d 300 -e $${IMG%%.svg}.png $$IMG; \
	done
	# render with pandoc
	pandoc \
		-o "mparker_thesis_draft_`date +%d%m%y`.pdf" \
		--template "template.tex" \
		-V fontsize=12pt \
		-V papersize=a4paper \
		-V documentclass=report \
		--pdf-engine=xelatex \
		introduction/tmp.md \
		chapter_3/tmp.md \
		chapter_4/tmp.md \
		chapter_5/tmp.md \
		chapter_6/tmp.md \
		discussion/tmp.md ;
	rm */tmp.md
	rm */figures/*.png
