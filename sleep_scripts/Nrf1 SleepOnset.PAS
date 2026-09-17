PROGRAM Delta;
{determine sleep onset after sleep deprivation and compare to visually assessed sleep onset
in a 2h window starting 0.5h before sleep deprivation end}

USES
  strings;
CONST
  h48   = 43200;
TYPE
  epoch = RECORD
            state          : char;
            bin            : array[0..400] of single;
            EEGv,EMGv,temp : single;
          END;
VAR
  DIRf,OUTf                : text;
  EEGf                     : file of epoch;
  spec4                    : epoch;
  str                      : string[2];
  fln5                     : string[5];
  fln                      : string;
  vs                       : char;
  t,m,hz,strain,numM,
  n,x,so_y,so              : integer;
  stop                     : boolean;
  star                     : array[1..1800] of integer;


BEGIN
  assign(DIRf,'Nrf1 Directory.txt');    reset(DIRf);
  assign(OUTf,'Nrf1_SleepOnset.out'); rewrite(OUTf);

  FOR strain:=1 TO 3 DO
  BEGIN
    readln(DIRf,str); readln(DIRf,numM);
    FOR m:=1 TO numM DO
    BEGIN
      readln(DIRf,fln5,so_y);
      {read visual assessed sleep onset; note that in the provided directory files the
      visually-assesed sleep onsets have been already replaced by the final consesus sleep onsets}
      fln:=trimright(fln5); writeln(str:3,m:3,fln:7);
      assign(EEGf,fln+str+'.smo'); reset(EEGf);

      seek(EEGf,h48+4950);
      FOR t:=1 TO 1800 DO         {read 2h worth of sleep-wake states starting 30 min before the end of the sleep deprivation}
      BEGIN
        read(EEGf,spec4); vs:=spec4.state;
        IF vs in ['n','2'] THEN star[t]:=1 ELSE star[t]:=0;
      END;
      close(EEGf);

      t:=0; n:=1; x:=0; stop:=false;
      REPEAT
        inc(t);
        IF (star[t]>0) AND (star[t+1]>0) THEN {if 2 consective 4s epoch are scored as NREM sleep start checking for a consolidated sleep bout}
        BEGIN
          n:=1; so:=t; stop:=false; x:=0;
          REPEAT
            inc(t); inc(n);
            IF star[t]<1 THEN inc(x);                         {counting interruptions by single 4s epochs; 2 in the first 22 4s epochs disqualify this bout as being 'consolidated'}
            IF (star[t]<1) AND (star[t-1]<1) THEN stop:=true; {2 consective 4s epochs not scored as NREM sleep in the first 22 4s epochs disqualify this bout as being 'consolidated'}
          UNTIL stop OR (x>2) OR (n>22) OR (t>=1800);
        END;
      UNTIL (n>22) OR (t>=1800);

      writeln(OUTf,str:4,m:4,so_y:6,(so+4950):6,so:6,(so_y-so-4950):6);
      {aligning visually scored sleep onset (i.e., 'so_y') counted from light onset on day 3 and the automatic sleep onset counted from 30 min before sleep deprivation end}
    END; {mouse}
  END;  {strain}
  close(DIRf);
  close(OUTf);
END.