# To setup, EITHER activate the bundled virtual environment tnxpersontime-venv, which has tnxpersontime and dependencies pre-installed:

# to activate the virtual environment tnxpersontime-venv
# make sure to cd into the tutorials directory, then:

source tnxpersontime-venv/bin/activate


# OR pip install the tnxpersontime python module (included in the root dir of this package) to your environment of choice

# activate the env of your choice
# cd into the tnxpersontime dir within the TNXpersontime package, then:

pip install .



### OPTIONAL FOR JUPYTER NOTEBOOK ###

# if running jupyter notebook and you to work within your venv - note that this creates a kernel on your machine:
# while the tnxpersontime-venv is activated:

ipython kernel install --user --name=tnxpersontime

# then in the jupyter notebook you can navigate to Kernel -> Change kernel -> tnxpersontime

# if you want to remove this kernel later:

ipython kernelspec remove tnxpersontime