module mod_types
  implicit none

  integer, parameter :: int8 = selected_int_kind(1)
  integer, parameter :: int16 = selected_int_kind(2)
  integer, parameter :: int32 = selected_int_kind(4)
  integer, parameter :: int64 = selected_int_kind(8)

  integer, parameter :: real32 = selected_real_kind(6, 37)
  integer, parameter :: real64 = selected_real_kind(15, 307)
end module mod_types
