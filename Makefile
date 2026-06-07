MODE = None
NAME = v01
all:
ifneq ($(NAME), v01)
	@(find . -type f -name "*" -print0 | xargs -0 sed -i'' -e "s/v01/$(NAME)/g")
	@(mv v01/v01.tex v01/$(NAME).tex)
	@(mv v01 $(NAME))
endif
	$(MAKE) -C $(NAME) MODE=$(MODE)
	cp $(NAME)/build/tex/$(NAME).pdf $(NAME)_rosenbaum_hikade.pdf

plots:
	$(MAKE) -C $(NAME) plot

clean:
	$(MAKE) -C $(NAME) clean

.PHONY: all clean
