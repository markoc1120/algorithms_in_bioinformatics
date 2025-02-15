from setuptools import find_packages, setup

setup(
    name="alignment_plus",
    version="0.1",
    packages=find_packages(),
    entry_points={
        "console_scripts": [
            "alignment_plus=main:main",
        ],
    },
    install_requires=[
        "numpy",
    ],
)
