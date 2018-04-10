################################################################################
#
#  Trace
#
################################################################################

doc"""
***
    trace(a::NfOrdElem) -> fmpz

> Returns the trace of $a$.
"""
function trace(a::AbsNfOrdAbsElem)
  return FlintZZ(trace(a.elem_in_nf))
end

