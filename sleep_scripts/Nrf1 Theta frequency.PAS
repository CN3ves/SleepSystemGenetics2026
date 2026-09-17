PROGRAM TPF;
{Theta-peak frequency TPF (and theta-peak power) in the 12h dark periods}
{Because TDW is prevalent in the dark and TPF is faster than in the light periods, }
{only the dark periods were use for determining TPF in TDW}
{Similarly, for REM sleep TPF, only light periods are considered}

USES
  strings;
CONST
  d1   = 21600;
  d4   = 86400;
TYPE
  epoch = RECORD
            state          : char;
            bin            : array[0..400] of single;
            EEGv,EMGv,temp : single;
          END;
VAR
  DIRf,OUTf,TDWf         : text;
  EEGf                   : file of epoch;
  spec4                  : epoch;
  str                    : string[2];
  fln5                   : string[5];
  t,m,hz,strain,so,day,
  numM,state,start,c     : integer;
  thetapeak,thetapow,
  localpeak,localpow,ref : single;
  vs                     : array[0..(d4+1)] of char;


BEGIN
  assign(DIRf,'Nrf1 Direcotory.txt'); reset(DIRf);
  assign(OUTf,'Nrf TDW TPF.out'); rewrite(OUTf);

  FOR strain:=1 TO 3 DO
  BEGIN
    readln(DIRf,str); readln(DIRf,numM);
    FOR m:=1 TO numM DO
    BEGIN
      readln(DIRf,fln5,so,ref); {use same reference as for sleep-wake EEG spectra states in baseline to normalize theta power}
	  fln:=trimright(fln5);
      
      {READ STATES}
      assign(TDWf,fln+'.tdw'); reset(TDWf); {files generated with Nrf1 Determine TDW.pas}
      FOR t:=1 TO d4 DO readln(TDWf,vs[t]);
      close(TDWf);
      vs[0]:=vs[1];
      vs[(d4+1)]:=vs[d4];

      {SET TO ZERO}
      c:=0; thetapow:=0.0; thetapeak:=0.0;

      {READ SPECTRA}
      assign(EEGf,fln+'.smo'); reset(EEGf);
      FOR day:=1 TO 4 DO
      BEGIN
        {start:=(day-1)*d1;}  {12h light if TPF in REM sleep}
        start:=10800+(day-1)*d1; {12h dark for TPF in TDW}
        seek(EEGf,start);
        FOR t:=start+1 TO start+10800 DO
        BEGIN
          read(EEGf,spec4);
          IF (vs[t-1] in ['w','9']) AND (vs[t]='9') AND (vs[t+1] in ['w','9']) THEN      {TDW}
          {IF (vs[t-1]='r') AND (vs[t]='r') AND (vs[t+1]='r') THEN}
          BEGIN
            inc(c);
            localpow:=0.0;
            FOR hz:=20 TO 50 DO {TDW} {16 TO 40 for REM sleep as theta is slower}
            BEGIN
              IF spec4.bin[hz]>localpow THEN
              BEGIN
                localpow:=spec4.bin[hz];
                localpeak:=hz;
              END;
            END;
            thetapow:=thetapow+localpow;
            thetapeak:=thetapeak+localpeak;
          END;
        END;
      END;
      close(EEGf);

      {WRITE OUTPUT}
      thetapow:=thetapow*1000000/c; thetapow:=thetapow*100/ref;
      thetapeak:=thetapeak/c;       thetapeak:=thetapeak/4;
      writeln(OUTf,str:3,m:4,thetapow:10:4,thetapeak:10:4,c:10);
      writeln(str:3,m:4);
    
	END; {mouse}
  END; {strain}
  
  close(DIRf);
  close(OUTf);
END.