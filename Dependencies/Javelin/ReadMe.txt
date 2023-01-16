 Javelin module (https://github.com/rosswhitfield/javelin) 
 extended with custom modules for input/output, analysis, 
 and custom interactions... 
 All src files for custom modules are located in "custom" 
 directory. 
----------------------------------------------------------
 Installation: 
 Activate conda env (e.g. "javelin")
 conda activate javelin 

 Install dependencies: 
 conda env update --file environment.yml

 Build using setuptools: 
 python setup.py build_ext --inplace 
 
 Install using pip: 
 pip install -e .


