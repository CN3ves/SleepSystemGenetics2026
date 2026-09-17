PROGRAM Process S;
{see 10.1073/pnas.0602006103 & 10.1523/JNEUROSCI.21-08-02610.2001}

USES
  strings;
CONST
  day4=86400;
  day1=21600;
TYPE
  epoch = Record
            state          : char;
            bin            : array[0..400] of single;
            EEGv,EMGv,temp : single;
          end;
VAR
  DIRf,OUTf                : text;
  SMOf                     : file of epoch;
  epoch4                   : epoch;
  str                      : string[2];
  fln5                     : string[5];
  fln                      : string;
  s                        : char;
  vs                       : array[1..day4] of byte;
  tim,local,refS,assp,assm,
  pS,pS0,pS1,dt,ti,td,
  facti,factd              : single;
  t,m,strain,numM,nc,so    : integer;


BEGIN
  pS0:=148.0;
  assp:=282.0; assm:=55.0; Ti:=7.9; Td:=1.9;
  dt:=4/3600; facti:=EXP(-dt/Ti); factd:=EXP(-dt/Td);
  {parameters as in 10.1073/pnas.0602006103
  upper, lower assymptotes in % of reference, time constants for increase and decrease in hours}

  assign(DIRf,'Nrf1_Directory.txt'); reset(DIRf);
  assign(OUTf,'Nrf1_ProcessS.out'); rewrite(OUTf);

  FOR strain:=1 TO 3 DO
  BEGIN
    readln(DIRf,str); readln(DIRf,numM);
    FOR m:=1 TO numM DO
    BEGIN
      readln(DIRf,fln5,so); fln:=trimright(fln5);
      assign(SMOf,fln+str+'.smo'); reset(SMOf);
      writeln(strain:3,m:4,fln+str:10);

	  {Read sleep-wake states for NREM sleep in entire 4day experiment}
      FOR t:=1 TO day4 DO
      BEGIN
        read(SMOf,epoch4);
        IF epoch4.state in ['n','2','5'] THEN vs[t]:=1 ELSE vs[t]:=2;
      END;
      close(SMOf);

     {set pS at t=0 to the average of the levels reached at the end of the 2 baselines, as the end does not depend on the intial value and steady state is assumed}
      pS:=pS0;
      FOR t:=1 TO (day1*2) DO
      BEGIN
        IF vs[t]=1 THEN pS:=assm+(pS-assm)*factd ELSE pS:=assp-(assp-pS)*facti; {Process S decreases in NREM sleep and increases in Waking and REM sleep}
        IF t=day1 THEN pS1:=pS;
      END;
      pS0:=(pS1+pS)/2;

     {determine 100% baseline reference for Process S as done for delta power; i.e., ZT8-12}
      pS:=pS0; refS:=0; nc:=0;
      FOR t:=1 TO (day1*2) DO
      BEGIN
        IF vs[t]=1 THEN pS:=assm+(pS-assm)*factd ELSE pS:=assp-(assp-pS)*facti;
        IF ((t>7200) AND (t<11000)) OR ((t>7200+day1) AND (t<11000+day1)) THEN
        BEGIN
          inc(nc);
          refS:=refS+pS;
        END;
      END;
      refS:=refS/nc;

      {simulate 4 days and write to output every 15 min}
      pS:=pS0; t:=0; tim:=0.0; local:=pS*100/refS;
      writeln(OUTf,str:4,m:4,tim:10:3,local:10:4,refS:10:4);
      FOR t:=1 TO day4 DO
      BEGIN
        IF vs[t]=1 THEN pS:=assm+(pS-assm)*factd ELSE pS:=assp-(assp-pS)*facti;
        IF t MOD 225 = 0 THEN
        BEGIN
          tim:=t; tim:=tim/900; local:=pS*100/refS;
          writeln(OUTf,str:4,m:4,tim:10:3,local:10:4,refS:10:4);
        END;
      END;
    END; {mouse}
  END; {strain}

  close(DIRf);
  close(OUTf);
END.