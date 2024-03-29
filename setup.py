from setuptools import setup, find_packages

version = '0.6.13'

setup(
    name="contactplot",
    version=version,
    author="Carlos Rios V.",
    author_email="crosvera@gmail.com",
    description="Contactplot parse DrSasa files and plots the atom/residue buried surfice",
    packages=find_packages(),
    py_modules=[],
    include_package_data=True,
    install_requires=[
        'matplotlib>=2.1',
        'numpy',
        'pandas',
        'docopt',
        'seaborn'
    ],
    zip_safe=True,
    scripts=['contactplot.py']
)

