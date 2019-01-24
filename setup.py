from setuptools import setup

setup(
    name='pySNPCall',
    url='https://github.com/baribefe/pySNPCallingPipeline/',
    author='Emmanuel Naziga',
    author_email='baribefe@gmail.com',
    packages=['pysnpcall'],
    install_requires=['numpy','pandas'],
    version='0.1',
    license='MIT',
    description='Python implementation of Java SNP calling pipeline (https://github.com/DSGlab/SNPCallingPipeline/)',
    # We will also need a readme eventually (there will be a warning)
    long_description=open('README.txt').read(),
)

