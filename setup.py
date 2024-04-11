from setuptools import setup, find_packages
import pathlib

here = pathlib.Path(__file__).parent.resolve()

# long_description = (here / "README.md").read_text(encoding="utf-8")

setup(
    name="starlight_toolkit",
    version="0.1.0",
    description="A set of tools to manage the Starlight spectral synthesis code.",
    # long_description=long_description,
    url="https://github.com/arielwrl/starlight_toolkit",
    author="Ariel Werle",
    author_email="ariel.werle@inaf.it",
    package_dir={"": "starlight_toolkit"},
    packages=find_packages(where="starlight_toolkit"),
    python_requires=">=3, <4",
    install_requires=['astropy>=5.2.2',
                      'matplotlib>=3.7.2',
                      'numpy>=1.17.4',
                      'scipy>=1.3.3',
                      'setuptools>=45.2.0'], 
    include_package_data=True
)

"""
python3 setup.py sdist bdist_wheel
twine upload dist/* --verbose
"""