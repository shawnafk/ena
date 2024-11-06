from setuptools import setup, find_packages

setup(
    name='gena',
    version='0.1.0',
    author='shawn',
    author_email='zhengjs@mail.ustc.edu.cn',
    description='GENA data processing package',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    url='https://github.com/jiangshanzheng/gena',
    packages=find_packages(),
    install_requires=open('requirements.txt').read().splitlines(),
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
    ],
    python_requires='>=3.7',
)
