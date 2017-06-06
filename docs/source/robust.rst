Robust Statistics
==========

.. function:: trim(x; prop=0.0, count=0)

  Return a copy of ``x`` with either ``count`` or proportion ``prop`` of the highest
  and lowest elements removed.  To compute the trimmed mean of ``x`` use
  ``mean(trim(x))``; to compute the variance use ``trimvar(x)``.

.. function:: winsor(x; prop=0.0, count=0)

  Return a copy of ``x`` with either ``count`` or proportion ``prop`` of the lowest
  elements of ``x`` replaced with the next-lowest, and an equal number of the
  highest elements replaced with the previous-highest.  To compute the Winsorized
  mean of ``x`` use ``mean(winsor(x))``.

.. function:: trimvar(x; prop=0.0, count=0)

  Compute the variance of the trimmed mean of ``x``. This function uses
  the Winsorized variance, as described in Wilcox (2010).
