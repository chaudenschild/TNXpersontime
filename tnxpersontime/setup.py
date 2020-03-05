from setuptools import setup

setup(name='tnxpersontime',
      version='0.1',
      description='Tool suite for risk and incidence calculations using TriNetX data',
      url='http://github.com/chaudenschild/TNXpersontime',
      author='Christian Haudenschild',
      author_email='christian.haudenschild@trinetx.com',
      license='MIT',
      py_modules=['tnxpersontime'],
      install_requires=['numpy', 'scipy', 'pandas'],
      zip_safe=False)
