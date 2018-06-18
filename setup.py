from setuptools import setup, find_packages

version = '0.6'

setup(
    name="contactplot",
    version=version,
    author="Carlos Rios V.",
    author_email="crios@cadetech.cl",
    description="DMBroker: DamageMonitor Message Broker",
    packages=find_packages(),
    py_modules=[],
    include_package_data=True,
    install_requires=[
        'matplotlib>=2.0',
        'numpy',
        'pandas',
        'docopt',
        'seaborn'
    ],
    zip_safe=True,
    scripts=['contactplot.py']
)

