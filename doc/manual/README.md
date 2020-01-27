# doc_template

LaTeX documentation template based on REVTeX (since we're aspiring for PRL
articles).

# Prerequisite Packages

To resolve dependencies for compilation, you'll need the following packages on 
top of a vanilla Ubuntu installation:

```bash
sudo apt-get install texlive texlive-science texlive-publishers
```

Once these have been installed, you should be able to simply issue `make`, and
you'll have a shiny pdf for your hard work.

# Standards

As old-fashioned as it may be, all files should have at most 80 characters per 
line. While most IDEs automatically line break at an appropriate position,
looking at raw files is enormously painful without adopting this convention.
