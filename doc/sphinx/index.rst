.. mash documentation master file, created by
   sphinx-quickstart on Fri May  8 14:15:00 2015.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

|

.. image:: rabbit.png

|

An efficient and versatile sketching library for biological sequences

|

Publication
============
RabbitSketch paper is under review.

.. .. toctree::
   :maxdepth: 1
   
..   data

Download
=========

* `Source Code <https://github.com/ZekunYin/RabbitSketch>`_

Install Python Package
=========================

Install from pypi:

.. code::

  pip3 install rabbitsketch --user

Install from source code:

.. code::

  git clone --recursive https://github.com/ZekunYin/RabbitSketch.git
  cd RabbitSketch
  pip3 install . --user

Build C++ Dynamic and Static Libraries
=======

.. code::

  git clone --recursive https://github.com/ZekunYin/RabbitSketch.git
  cd RabbitSketch
  mkdir build && cd build
  cmake -DBUILDCXX=ON ..
  make -j4

Documentation
=============

.. toctree::
   :maxdepth: 2
   
   tutorials
   sketches
   distances

