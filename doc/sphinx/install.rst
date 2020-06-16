Install 
=======

Download
---------

* `Source Code <https://github.com/ZekunYin/RabbitSketch>`_

Install Python Package
-------------------------

Install from pypi:

.. code::

  pip3 install rabbitsketch --user

Install from source code:

.. code::

  git clone --recursive https://github.com/ZekunYin/RabbitSketch.git
  cd RabbitSketch
  pip3 install . --user

Build C++ Dynamic and Static Libraries
----------------------------------------

.. code::

  git clone --recursive https://github.com/ZekunYin/RabbitSketch.git
  cd RabbitSketch
  mkdir build && cd build
  cmake -DBUILDCXX=ON ..
  make -j4
  make install //todo change it to make install 
