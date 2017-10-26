.. highlight:: rst

Shape
=====
This is the shape module

.. toctree::
  :hidden:

  shape_references

.. automodule:: symeess.shape
   :members:

.. py:class:: symeess.shape

   .. py:method:: get_measure(geometry)
      Compute the shape measure of the given geometry

      :param obj geometry: this object contains information about the positions of the molecule as well as the
                           information needed to carry out the shape's module
      :return: difference between user's structure and the one's that is compared. 0 is no difference and 100
               is completly different
   .. py:method:: get_structure(geometry)
      Calculate the ideal structure of the given geometry from a reference structure

      :param obj geometry: same as before
      :return: ideal structure if user's structure had the compared structure's shape

:doc:`shape_references`
