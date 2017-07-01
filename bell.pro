;+
; NAME:
;       Bell
;
; PURPOSE:
;
;       This program plays a tone (or any *.wav file) on a Windows machine.
;
; AUTHOR:
;
;       FANNING SOFTWARE CONSULTING
;       David Fanning, Ph.D.
;       1645 Sheely Drive
;       Fort Collins, CO 80526 USA
;       Phone: 970-221-0438
;       E-mail: davidf@dfanning.com
;       Coyote's Guide to IDL Programming: http://www.dfanning.com
;
; CATEGORY:
;
;       General Programming
;
; CALLING SEQUENCE:
;
;       Bell
;
; OPTIONAL INPUTS:
;
;       None.
;
; OUTPUTS:
;
;       None. A tone or sound is played.
;
; INPUT KEYWORDS:
;
;        ALLWRONG  -- Set this keyword to play the file "allwrong.wav". (Comes in the zip file.)
;        CHIMES    -- Set this keyword to play the file "chimes.wav".
;        CHORD     -- Set this keyword to play the file "chord.wav".
;        COYOTE    -- Set this keyword to play the file "coyote.wav".
;        DING      -- Set this keyword to play the file "ding.wav". This is the default.
;        DOH       -- Set this keyword to play the file "doh.wav". (Comes in the zip file.)
;        FILENAME  -- Set this keyword to the name of a *.wav file to play.
;        NOTIFY    -- Set this keyword to play the file "notify.wav".
;        TADA      -- Set this keyword to play the file "tada.wav".
;        WARNING   -- Set this keyword to play the file "warning.wav". (Comes in the zip file.)
;
;        Any named *.wav file that doesn't come with this distribution will be found in
;        the Windows Media directory. This is normally C:\WINNT\MEDIA.
;
; OUTPUT KEYWORDS:
;
;       None.
;
; DEPENDENCIES:
;
;       This is a Windows-only application and files have to be located in
;       the correct place for them to work.
;
;       The file winclip.dll has to go into the main IDL directory for *each*
;       of your IDL distributions. That is "C:\RSI\IDL5.5" or "C:\RSI\IDL5.6"
;       or wherever you have installed IDL.
;
;       You can put the *.wav files anywhere on your IDL path or in the
;       main Windows Media directory (usually "C:\WINNT\MEDIA").
;
;       Put the bell.pro file anywhere in your IDL path.
;
; MODIFICATION HISTORY:
;
;       Thrown together by David W. Fanning from files supplied by Peter Mason
;       and Andrew Cool, two Aussies who know how to have fun! 14 Nov 2002.
;-
;
;###########################################################################
;
; LICENSE
;
; This software is OSI Certified Open Source Software.
; OSI Certified is a certification mark of the Open Source Initiative.
;
; Copyright © 2002 Fanning Software Consulting.
;
; This software is provided "as-is", without any express or
; implied warranty. In no event will the authors be held liable
; for any damages arising from the use of this software.
;
; Permission is granted to anyone to use this software for any
; purpose, including commercial applications, and to alter it and
; redistribute it freely, subject to the following restrictions:
;
; 1. The origin of this software must not be misrepresented; you must
;    not claim you wrote the original software. If you use this software
;    in a product, an acknowledgment in the product documentation
;    would be appreciated, but is not required.
;
; 2. Altered source versions must be plainly marked as such, and must
;    not be misrepresented as being the original software.
;
; 3. This notice may not be removed or altered from any source distribution.
;
; For more information on Open Source Software, visit the Open Source
; web site: http://www.opensource.org.
;
;###########################################################################


;============================================================
; Routines for Windows stuff.
; Peter Mason, CSIRO D.E.M., 1996
; e-mail: peter.mason@dem.csiro.au
; These routines use the DLL:  WINCLIP.DLL which must be in the IDL installation dir
;============================================================

function wccheck,shlib
; Check that we'll be able to use the procedures herein
; The function return value is 0 if OK, 1 if not.
; Also, the path\filename of the WINCLIP DLL is returned in SHLIB if all's ok
; To make life simple, I insist on having winclip.dll stored in the main IDL directory.
; It doesn't really have to be, of course, but if you put it somewhere else then you will
; have to modify this routine accordingly.

shlib=''
if (strpos(strupcase(!version.os),'WIN') ne 0) then begin
  message,/cont,'Copying to the clipboard is only supported under Windows platforms.'
  return,1
endif
shlib=findfile(!dir+'\winclip.dll')
shlib=strtrim(shlib(0),2)
if (shlib eq '') then begin
  message,/cont,'Copying to the clipboard requires '+!dir+'\WINCLIP.DLL'
  return,1
endif
return,0
end

;============================================================
; This plays a sound (asynchronously) from a .WAV file.
; The first parameter is a scalar string specifying the .WAV filename.
; The second, optional parameter gives a couple of options:
;   0    Play the sound once (already!) - the default
;   1    Play the sound repeatedly - in a loop
;   2    Stop playing sounds for IDL (good with option 1 :-)
pro playwavfile,wavfile,opt=opt
if wccheck(shlib) then return
b=[byte(wavfile),0b]          ;sorry I just can't remember the convolutions with call_ex and strings - a null-terminated bytarr is simpler!
if n_elements(opt) eq 1 then opt=(long(opt)>0L)<2L $
else opt=0L
j=call_external(shlib,'idlplaywav',b,opt)
end


FUNCTION BELL_ERROR_MESSAGE, theMessage, Error=error, Informational=information, $
   Traceback=traceback, NoName=noname, Title=title, _Extra=extra

On_Error, 2

   ; Check for presence and type of message.

IF N_Elements(theMessage) EQ 0 THEN theMessage = !Error_State.Msg
s = Size(theMessage)
messageType = s[s[0]+1]
IF messageType NE 7 THEN BEGIN
   Message, "The message parameter must be a string.", _Extra=extra
ENDIF

   ; Get the call stack and the calling routine's name.

Help, Calls=callStack
IF Float(!Version.Release) GE 5.2 THEN $
   callingRoutine = (StrSplit(StrCompress(callStack[1])," ", /Extract))[0] ELSE $
   callingRoutine = (Str_Sep(StrCompress(callStack[1])," "))[0]

   ; Are widgets supported?

widgetsSupported = ((!D.Flags AND 65536L) NE 0)
IF widgetsSupported THEN BEGIN

      ; If this is an error produced with the MESSAGE command, it is a trapped
      ; error and will have the name "IDL_M_USER_ERR".

   IF !ERROR_STATE.NAME EQ "IDL_M_USER_ERR" THEN BEGIN

      IF N_Elements(title) EQ 0 THEN title = 'Trapped Error'

         ; If the message has the name of the calling routine in it,
         ; it should be stripped out. Can you find a colon in the string?

      colon = StrPos(theMessage, ":")
      IF colon NE -1 THEN BEGIN

            ; Extract the text up to the colon. Is this the same as
            ; the callingRoutine? If so, strip it.

         IF StrMid(theMessage, 0, colon) EQ callingRoutine THEN $
            theMessage = StrMid(theMessage, colon+1)

      ENDIF

         ; Add the calling routine's name, unless NONAME is set.

      IF Keyword_Set(noname) THEN BEGIN
         answer = Dialog_Message(theMessage, Title=title, _Extra=extra, $
            Error=error, Information=information)
      ENDIF ELSE BEGIN
         answer = Dialog_Message(StrUpCase(callingRoutine) + ": " + $
            theMessage, Title=title, _Extra=extra, $
            Error=error, Information=information)
      ENDELSE

   ENDIF ELSE BEGIN

         ; Otherwise, this is an IDL system error.

      IF N_Elements(title) EQ 0 THEN title = 'System Error'

      IF StrUpCase(callingRoutine) EQ "$MAIN$" THEN $
         answer = Dialog_Message(theMessage, _Extra=extra, Title=title, $
            Error=error, Information=information) ELSE $
      IF Keyword_Set(noname) THEN BEGIN
         answer = Dialog_Message(theMessage, _Extra=extra, Title=title, $
            Error=error, Information=information)
      ENDIF ELSE BEGIN
         answer = Dialog_Message(StrUpCase(callingRoutine) + "--> " + $
            theMessage, _Extra=extra, Title=title, $
            Error=error, Information=information)
      ENDELSE
   ENDELSE
ENDIF ELSE BEGIN
      Message, theMessage, /Continue, /NoPrint, /NoName, /NoPrefix, _Extra=extra
      Print, '%' + callingRoutine + ': ' + theMessage
      answer = 'OK'
ENDELSE

   ; Provide traceback information if requested.

IF Keyword_Set(traceback) THEN BEGIN
   Help, /Last_Message, Output=traceback
   Print,''
   Print, 'Traceback Report from ' + StrUpCase(callingRoutine) + ':'
   Print, ''
   FOR j=0,N_Elements(traceback)-1 DO Print, "     " + traceback[j]
ENDIF

RETURN, answer
END ; ----------------------------------------------------------------------------


PRO Bell, $
   AllWrong=allwrong, $ ; allwrong.wav
   Chimes=chimes, $     ; chimes.wav
   Chord=chord, $       ; chord.wav
   Coyote=coyote, $     ; coyote.wav
   Ding=ding, $         ; ding.wav -- The default.
   Doh=doh, $           ; doh.wav
   Filename=filename, $ ; The name of a "wav" file to play.
   Notify=notify, $     ; notify.wav
   Tada=tada, $         ; tada.wav
   Warning=warning      ; warning.wav

Catch, theError
IF theError NE 0 THEN BEGIN
   Catch, /Cancel
   ok = Bell_Error_Message(/Traceback)
   RETURN
ENDIF

theToneFile = 'ding.wav'
IF Keyword_Set(ding) THEN theToneFile = 'ding.wav'
IF Keyword_Set(chord) THEN theToneFile = 'chord.wav'
IF Keyword_Set(chimes) THEN theToneFile = 'chimes.wav'
IF Keyword_Set(coyote) THEN theToneFile = 'coyote.wav'
IF Keyword_Set(doh) THEN theToneFile = 'doh.wav'
IF Keyword_Set(notify) THEN theToneFile = 'notify.wav'
IF Keyword_Set(tada) THEN theToneFile = 'tada.wav'
IF Keyword_Set(warning) THEN theToneFile = 'warning.wav'
IF Keyword_Set(allwrong) THEN theToneFile = 'allwrong.wav'
IF N_Elements(filename) NE 0 THEN theToneFile = filename

   ; Can you find the file?

file = File_Which(theToneFile, /INCLUDE_CURRENT_DIR)
IF file EQ "" THEN file = File_Which('C:\WINNT\MEDIA\',theToneFile)
IF file EQ "" THEN BEGIN
   Message, 'Cannot find the "' + theToneFile +'" file.', /Informational
   BEEP
ENDIF ELSE Playwavfile, file, OPT = 0

END