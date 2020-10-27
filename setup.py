from setuptools import find_packages, setup

requirements = [
    "azcam",
    "numpy",
    "scipy",
    "PyPDF2",
    "pdfkit",
    "markdown",
]


with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name="azcam-testers",
    version="20.1",
    description="azcam extension for sensor characterization",
    long_description=long_description,
    author="Michael Lesser",
    author_email="mlesser@arizona.edu",
    keywords="ccd imaging astronomy observation observatory",
    packages=find_packages(),
    zip_safe=False,
    install_requires=requirements,
)
