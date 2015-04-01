from setuptools import setup, find_packages
setup(
    name = "aospy",
    version = "0.0",
    packages = find_packages(),
    scripts = ['main.py', 'print_table.py'],

    author = "Spencer A. Hill",
    author_email = "spencerahill@gmail.com",
    description = "Automated gridded climate data analysis and visualization",
    license = "GPL",
    keywords = "climate science",
    url = "https://github.com/spencerahill/aospy"
)
