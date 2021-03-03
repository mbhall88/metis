"""The setup script."""

from setuptools import setup, find_packages

with open('README.md') as readme_file:
    readme = readme_file.read()

requirements = ['Click>=6.0', 'numpy', 'matplotlib', 'seaborn', 'six', 'pysam']

setup_requirements = ['pytest-runner', ]

test_requirements = ['pytest', 'click', 'pysam', 'six',
                     'matplotlib']

setup(
    author="Michael Benjamin Hall",
    author_email='mbhall88@gmail.com',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Programming Language :: Python :: 3.6',
    ],
    description="Quality control plotting for long reads",
    entry_points={
        'console_scripts': [
            'pistis=pistis.pistis:main',
        ],
    },
    install_requires=requirements,
    license="MIT license",
    long_description=readme,
    long_description_content_type="text/markdown",
    include_package_data=True,
    keywords='pistis',
    name='pistis',
    packages=find_packages(include=['pistis']),
    setup_requires=setup_requirements,
    test_suite='tests',
    tests_require=test_requirements,
    url='https://github.com/mbhall88/pistis',
    version='0.3.4',
    zip_safe=False,
)
