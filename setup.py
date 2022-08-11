from setuptools import setup, find_packages

setup(name='bihidef',
    version='0.1',
    description='Bipartite adaptation of HiDef using BRIM',
    url='https://github.com/genisott/BiHiDef',
    author='Genís Calderer',
    author_email='genis.calderer@gmail.com',
    license='MIT',
    packages=['bihidef'],
    install_requires=['pandas',
    'numpy',
    'python-igraph',"networkx","scipy","netZooPy","hidef"],
    zip_safe=False)
