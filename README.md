To run go-no-go master

conda create -n py37 python=3.7
conda activate py37
python -m venv env
env/Scripts/activate.bat
preinstall.bat
python -m pip install -e .



to install, and


conda activate py37
env\Scripts\activate.bat
python main.py
env\Scripts\deactivate.bat


to run (edited) 


deactivate:

env\Scripts\deactivate.bat


To install toon:

python -m pip install toon==0.12.8

To install pip:
python -m pip install pip==20.0