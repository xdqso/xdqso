;+
; NAME:
;   xdqso_version
; PURPOSE:
;   Return the version name for the xdqso product
; CALLING SEQUENCE:
;   vers = xdqso_version()
; OUTPUTS:
;   vers       - Version name for the product xdqso
; COMMENTS:
;   Requires that the XDQSO_DIR environment variable be set
; VERSION:
;   $Id$
;-
;------------------------------------------------------------------------------
FUNCTION xdqso_version
    RETURN, FILE_BASENAME(GETENV('XDQSO_DIR'))
END
;------------------------------------------------------------------------------
