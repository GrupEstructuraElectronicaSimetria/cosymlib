.. highlight:: rst

Shape
=====
This is the shape module

.. toctree::
  :hidden:

  shape_references

.. py:class:: symeess.shape

   .. py:method:: get_measure(geometry)
      Method that prints to file shape's measure

      :param obj geometry: his object contains information about the positions of the molecule as well as the information
                              needed to carry out the shape's module
      :return: difference between user's structure and the one's that is compared. While 0 is no difference and 100
               is completly different
   .. py:method:: write_shape_structure_2file(shape_label, central_atom=None, output_name='../examples/symeess_shape')
       Method that prints to file shape's structure

      :param obj geometry: same as before
      :return: ideal structure if user's structure had the compared structure's shape

:doc:`shape_references`
