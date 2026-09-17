PROGRAM SleepWake;
{hourly values of time spent in sleep-wake states over the 4day experiment}

USES
  strings;
TYPE
  epoch = RECORD
            state          : char;
            bin            : array[0..400] of single;
            EEGv,EMGv,temp : single;
          END;
VAR
  DIRf,OUTf                : text;
  SMOf                     : file of epoch;
  fln5                     : string[5];
  fln                      : string;
  str                      : string[2];
  m,numM,h,t,wc,nc,rc,
  strain,day,start         : integer;
  spec4                    : epoch;
  wr,nr,rr,wacc,nacc,racc,
  tim                      : single;
  w,n,r                    : array[1..4,1..24] of single;


BEGIN
  assign(DIRf,'Nrf1 Directory.txt');    reset(DIRf);
  assign(OUTf,'Nrf1_SWstates1h.out'); rewrite(OUTf);

  FOR strain:=1 TO 3 DO
  BEGIN
    readln(DIRf,str); readln(DIRf,numM);

    FOR m:=1 TO numM DO
    BEGIN
      readln(DIRf,fln5); fln:=trimright(fln5); writeln(str:3,m:3,fln:7);
      assign(SMOf,fln+str+'.smo'); reset(SMOf);

      FOR day:=1 TO 4 DO
      BEGIN
        FOR h:=1 TO 24 DO
        BEGIN
          wc:=0; nc:=0; rc:=0;
          FOR t:=1 TO 900 DO {900 4s epcochs = 1h}
          BEGIN
            read(SMOf,spec4);
            CASE spec4.state OF
              'w','1','4' : inc(wc);
              'n','2','5' : inc(nc);
              'r','3','6' : inc(rc);
            END;
          END; {t}
          wr:=wc; w[day,h]:=wr/15; {minutes per hour}
          nr:=nc; n[day,h]:=nr/15;
          rr:=rc; r[day,h]:=rr/15;
        END; {h}
      END; {day}
      close(SMOf);

      {average the two baseline days}
      FOR h:=1 TO 24 DO
      BEGIN
        w[2,h]:=(w[1,h]+w[2,h])/2;
        n[2,h]:=(n[1,h]+n[2,h])/2;
        r[2,h]:=(r[1,h]+r[2,h])/2;
      END;

      {write waking, NREM sleep, and REM sleep values to output}
      {day=2: average of 2 baseline, day=3: sleep deprivation+recovery day 1, day=4: recovery day 2}
      FOR day:=2 TO 4 DO
      BEGIN
        FOR h:=1 TO 24 DO
        BEGIN
          tim:=h; tim:=(tim-0.5)+((day-2)*24); {hourly values must be plotted at interval midpoint}
          writeln(OUTf,str:3,m:4,day:3,h:3,tim:5:1,w[day,h]:10:4,n[day,h]:10:4,r[day,h]:10:4);
        END;
      END

    END; {mouse}
  END; {strain}

  close(DIRf);
  close(OUTf);
END.