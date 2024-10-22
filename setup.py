import os
from setuptools import setup, find_packages


PACKAGE = "RenShape"

def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()
    
def find_scripts(path):
    temp = [path+'/'+s for s in os.listdir(path)] 
    try:
        temp.remove(path+'/'+'README.md')
    except:
        pass
    return temp


setup( # Finally, pass this all along to distutils to do the heavy lifting.
    name             = PACKAGE,
    version          = '1.0.0',
    description      = 'RenShape packages',
    long_description = read('README.md'),
    author           = ' ',
    author_email     = ' ',
    url              = 'https://github.com/matteo-brg/RenShape',
#    license          = read('LICENSE'),
    install_requires = read('requirements.txt').splitlines(),
    package_dir      = {'': 'src'},
    packages         = find_packages(where='src'),
    scripts          = find_scripts('scripts'),
#    include_package_data=True,
    python_requires  = '>=3.8', 
    zip_safe         = False,
)
