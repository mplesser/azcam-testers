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
    version="21.1.2",
    description="azcam extension for sensor characterization",
    long_description=long_description,
    long_description_content_type="text/markdown",
    author="Michael Lesser",
    author_email="mlesser@arizona.edu",
    keywords="ccd imaging astronomy observation observatory",
    packages=find_packages(),
    zip_safe=False,
    install_requires=requirements,
)
