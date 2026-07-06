c     ******************************************************************
c     * functions.h                                                    *
c     * contains all of the function type declarations in towhee       *
c     *                                                                *
c     * originally written 11-12-2002 by M.G. Martin                   *
c     * last modified 03-30-2020 by M.G. Martin                        *
c     ******************************************************************
#ifdef FUNCTION_ARCCOS
      double precision twh_arccos
#undef FUNCTION_ARCCOS
#endif
c
#ifdef FUNCTION_CHECK_COMMENT
      logical twh_check_comment
#undef FUNCTION_CHECK_COMMENT
#endif
c
#ifdef FUNCTION_CHECK_LABEL
      logical twh_check_label
#undef FUNCTION_CHECK_LABEL
#endif
c
#ifdef FUNCTION_CMP_EQ
#ifdef SAFE_COMPARE
      logical twh_cmp_eq
#else
#define twh_cmp_eq(A,B) ((A).eq.(B))
#endif
#undef FUNCTION_CMP_EQ
#endif
c
#ifdef FUNCTION_CMP_GT
#ifdef SAFE_COMPARE
      logical twh_cmp_gt
#else
#define twh_cmp_gt(A,B) ((A).gt.(B))
#endif
#undef FUNCTION_CMP_GT
#endif
c
#ifdef FUNCTION_CMP_LT
#ifdef SAFE_COMPARE
      logical twh_cmp_lt
#else
#define twh_cmp_lt(A,B) ((A).lt.(B))
#endif
#undef FUNCTION_CMP_LT
#endif
c
#ifdef FUNCTION_DERF
      double precision twh_derf 
#undef FUNCTION_DERF
#endif
c
#ifdef FUNCTION_DERFC
      double precision twh_derfc
#undef FUNCTION_DERFC
#endif
c
#ifdef FUNCTION_DISTANCE
      double precision twh_distance
#undef FUNCTION_DISTANCE
#endif
c
#ifdef FUNCTION_DOTPRODUCT
      double precision twh_dotproduct
#undef FUNCTION_DOTPRODUCT
#endif
c
#ifdef FUNCTION_EWALD_CORRECT
      double precision twh_ewald_correct
#undef FUNCTION_EWALD_CORRECT
#endif
c
#ifdef FUNCTION_EWALD_SELF
      double precision twh_ewald_self
#undef FUNCTION_EWALD_SELF
#endif
c
#ifdef FUNCTION_EXPON
      double precision twh_expon 
#undef FUNCTION_EXPON
#endif
c
#ifdef FUNCTION_EXTRACTDENS
      double precision twh_extractdens 
#undef FUNCTION_EXTRACTDENS
#endif
c
#ifdef FUNCTION_EXTRACTEMBED
      double precision twh_extractembed
#undef FUNCTION_EXTRACTEMBED
#endif
c
#ifdef FUNCTION_EXTRACTPAIR
      double precision twh_extractpair 
#undef FUNCTION_EXTRACTPAIR
#endif
c
#ifdef FUNCTION_FEBIAS
      double precision twh_febias 
#undef FUNCTION_FEBIAS
#endif
c
#ifdef FUNCTION_GAUSSPROB
      double precision twh_gaussprob 
#undef FUNCTION_GAUSSPROB
#endif
c
#ifdef FUNCTION_GAUSSIAN
      double precision twh_gaussian 
#undef FUNCTION_GAUSSIAN
#endif
c
#ifdef FUNCTION_GETATOMNUM
      integer twh_getatomnum
#undef FUNCTION_GETATOMNUM
#endif
c
#ifdef FUNCTION_GETNBTYPE
      integer twh_getnbtype
#undef FUNCTION_GETNBTYPE
#endif
c
#ifdef FUNCTION_GETSIGN
      double precision twh_getsign 
#undef FUNCTION_GETSIGN
#endif
c
#ifdef FUNCTION_GET_ATOMLIST_CORD
      integer twh_get_atomlist_cord
#undef FUNCTION_GET_ATOMLIST_CORD
#endif
c
#ifdef FUNCTION_GET_DERIVATIVE
      double precision scp_get_derivative
#undef FUNCTION_GET_DERIVATIVE
#endif
c
#ifdef FUNCTION_GET_FOREIGN_LAMBDA
      double precision scp_get_foreign_lambda
#undef FUNCTION_GET_FOREIGN_LAMBDA
#endif
c
#ifdef FUNCTION_GET_LAAHERE
      logical twh_get_laahere
#undef FUNCTION_GET_LAAHERE
#endif
c
#ifdef FUNCTION_GET_LBENDHERE
      logical twh_get_lbendhere
#undef FUNCTION_GET_LBENDHERE
#endif
c
#ifdef FUNCTION_GET_LHERE
      logical twh_get_lhere
#undef FUNCTION_GET_LHERE
#endif
c
#ifdef FUNCTION_GET_LIMPHERE
      logical twh_get_limphere
#undef FUNCTION_GET_LIMPHERE
#endif
c
#ifdef FUNCTION_GET_LOFHERE
      logical twh_get_lofhere
#undef FUNCTION_GET_LOFHERE
#endif
c
#ifdef FUNCTION_GET_LOFTOR
      logical twh_get_loftor
#undef FUNCTION_GET_LOFTOR
#endif
c
#ifdef FUNCTION_GET_LTORHERE
      logical twh_get_ltorhere
#undef FUNCTION_GET_LTORHERE
#endif
c
#ifdef FUNCTION_GET_LVIBHERE
      logical twh_get_lvibhere
#undef FUNCTION_GET_LVIBHERE
#endif
c
#ifdef FUNCTION_GET_NATIVE_LAMBDA
      double precision scp_get_native_lambda
#undef FUNCTION_GET_NATIVE_LAMBDA
#endif
c
#ifdef FUNCTION_GET_SCALING_STYLE
      integer scp_get_scaling_style
#undef FUNCTION_GET_SCALING_STYLE
#endif
c
#ifdef FUNCTION_GET_SCALING_STYLE_STRING
      character*30 scp_get_scaling_style_string
#undef FUNCTION_GET_SCALING_STYLE_STRING
c
#endif
c
#ifdef FUNCTION_GET_STRING_LENGTH
      integer twh_get_string_length
#undef FUNCTION_GET_STRING_LENGTH
#endif
c
#ifdef FUNCTION_IN_ATOMLIST
       logical twh_in_atomlist
#undef FUNCTION_IN_ATOMLIST
#endif
c
#ifdef FUNCTION_INVERSELAWOFCOSINE
      double precision twh_inverselawofcosine
#undef FUNCTION_INVERSELAWOFCOSINE
#endif
c
#ifdef FUNCTION_LAWOFCOSINE
      double precision twh_lawofcosine
#undef FUNCTION_LAWOFCOSINE
#endif
c
#ifdef FUNCTION_LEN_TRIM
      integer twh_len_trim
#undef FUNCTION_LEN_TRIM
#endif
c
#ifdef FUNCTION_LINCLUDE
      logical twh_linclude
#undef FUNCTION_LINCLUDE
#endif
c
#ifdef FUNCTION_MAXBOXLENGTH
      double precision twh_maxboxlength 
#undef FUNCTION_MAXBOXLENGTH
#endif
c
#ifdef FUNCTION_MINBOXLENGTH
      double precision twh_minboxlength
#undef FUNCTION_MINBOXLENGTH
#endif
c
#ifdef FUNCTION_OLDGETATOMNUM
      integer twh_oldgetatomnum
#undef FUNCTION_OLDGETATOMNUM
#endif
c
#ifdef FUNCTION_PEEK_LABEL
      character*50 twh_peek_label
#undef FUNCTION_PEEK_LABEL
#endif
c
#ifdef FUNCTION_RANDOM
      double precision twh_random 
#undef FUNCTION_RANDOM
#endif
c
#ifdef FUNCTION_READ_DIR_STRING
      character*MAXDIRLENGTH twh_read_dir_string
#undef FUNCTION_READ_DIR_STRING
#endif
c
#ifdef FUNCTION_READ_FLOAT
      double precision twh_read_float
#undef FUNCTION_READ_FLOAT
#endif
c
#ifdef FUNCTION_READ_INTEGER
      integer twh_read_integer
#undef FUNCTION_READ_INTEGER
#endif
c
#ifdef FUNCTION_READ_LABELED_FLOAT
      double precision twh_read_labeled_float
#undef FUNCTION_READ_LABELED_FLOAT
#endif
c
#ifdef FUNCTION_READ_LABELED_INTEGER
      integer twh_read_labeled_integer
#undef FUNCTION_READ_LABELED_INTEGER
#endif
c
#ifdef FUNCTION_READ_LABELED_LOGICAL
      logical twh_read_labeled_logical
#undef FUNCTION_READ_LABELED_LOGICAL
#endif
c
#ifdef FUNCTION_READ_LOGICAL
      logical twh_read_logical
#undef FUNCTION_READ_LOGICAL
#endif
c
#ifdef FUNCTION_READ_STRING
      character*(50) twh_read_string
#undef FUNCTION_READ_STRING
#endif
c
#ifdef FUNCTION_SCALE_ATOMS
      logical scp_scale_atoms
#undef FUNCTION_SCALE_ATOMS
#endif
c
#ifdef FUNCTION_LIMITED_DOUBLE
      double precision twh_limited_double
#undef FUNCTION_LIMITED_DOUBLE
#endif
c
#ifdef FUNCTION_SAFE_DOUBLE
      double precision twh_safe_double
#undef FUNCTION_SAFE_DOUBLE
#endif
c
#ifdef FUNCTION_VANGANG
      double precision twh_vangang 
#undef FUNCTION_VANGANG
#endif
c
#ifdef FUNCTION_VANGLE
      double precision twh_vangle 
#undef FUNCTION_VANGLE
#endif
c
#ifdef FUNCTION_VBOND
      double precision twh_vbond 
#undef FUNCTION_VBOND
#endif
c
#ifdef FUNCTION_VBONBON
      double precision twh_vbonbon 
#undef FUNCTION_VBONBON
#endif
c
#ifdef FUNCTION_VEEFONE
      double precision twh_veefone 
#undef FUNCTION_VEEFONE
#endif
c
#ifdef FUNCTION_VEMBED
      double precision twh_vembed
#undef FUNCTION_VEMBED
#endif
c
#ifdef FUNCTION_VFIELD
      double precision twh_vfield
#undef FUNCTION_VFIELD
#endif
c
#ifdef FUNCTION_VIMPROPER
      double precision twh_vimproper 
#undef FUNCTION_VIMPROPER
#endif
c
#ifdef FUNCTION_VONEFIVE
      double precision twh_vonefive
#undef FUNCTION_VONEFIVE
#endif
c
#ifdef FUNCTION_VSASA
      double precision twh_vsasa 
#undef FUNCTION_VSASA
#endif
c
#ifdef FUNCTION_VTHREEBODY
      double precision twh_vthreebody
#undef FUNCTION_VTHREEBODY
#endif
c
#ifdef FUNCTION_VTORSION
      double precision twh_vtorsion
#undef FUNCTION_VTORSION
#endif
c
#ifdef FUNCTION_VTWOBODY
      double precision twh_vtwobody 
#undef FUNCTION_VTWOBODY
#endif
c
#ifdef FUNCTION_WCOULOMB
      double precision twh_wcoulomb
#undef FUNCTION_WCOULOMB
#endif
c
#ifdef FUNCTION_WTWOBODY
      double precision twh_wtwobody
#undef FUNCTION_WTWOBODY
#endif
