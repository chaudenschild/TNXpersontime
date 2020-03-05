#create virtual env
python3 -m venv env

#activate venv
source env/bin/activate

#install requirements
pip install numpy pandas scipy

#install local version of tnxpersontime python module to venv
cd tnxpersontime
pip install .
