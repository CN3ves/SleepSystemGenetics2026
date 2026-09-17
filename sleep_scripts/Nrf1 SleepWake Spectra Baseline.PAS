PROGRAM EEG spectra;
{for each sleep-wake state in 48h baseline}

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
  DIRf,OUTf            : text;
  EEGf                 : file of epoch;
  spec4                : epoch;
  str                  : string[2];
  fln5                 : string[5];
  fln                  : string;
  t,m,hz,strain,
  numM,state           : integer;
  c                    : array[1..3] of integer;
  hzr,local,reft       : single;
  rat,ref              : array[1..3] of single;
  spec_state           : array[1..3,3..360] of single;
  star                 : array[0..(h48+1)] of char;


BEGIN
  assign(DIRf,'Nrf1_Directory.txt');    reset(DIRf);
  assign(OUTf,'Nrf1_SW_Spectra.out'); rewrite(OUTf);

  {Wake, NREM sleep, REM sleep fractions in baseline in 2 baseline of the entire experiment, i.e. over 31 mice, calculated elsewhere}
  {to weigh the reference according to the prevalence of each state independent of eventual group differences in sleep-wake times}
  rat[1]:=0.52743; rat[2]:=0.42195; rat[3]:=0.05098;

  FOR strain:=1 TO 3 DO
  BEGIN
    readln(DIRf,str); readln(DIRf,numM);
    FOR m:=1 TO numM DO
    BEGIN
      {set spectral arrays to zero}
      FOR state:=1 TO 3 DO
      BEGIN
        c[state]:=0;
        FOR hz:=3 TO 360 DO spec_state[state,hz]:=0.0;
      END;

      readln(DIRf,fln5); fln:=trimright(fln5); writeln(str:3,m:3,fln:7);
      assign(EEGf,fln+str+'.smo'); reset(EEGf);

      {read states in 2 baseline days}
      FOR t:=1 TO h48 DO
      BEGIN
        read(EEGf,spec4);
        star[t]:=spec4.state;
      END;
      star[0]:=star[1];
      star[(h48+1)]:=star[h48];

      {read spectra}
      reset(EEGf);
      FOR t:=1 TO h48 DO
      BEGIN
        read(EEGf,spec4);
        IF spec4.state in ['w','n','r'] THEN
        BEGIN
          IF (star[t-1]=spec4.state) AND (star[t+1]=spec4.state) THEN 
		  {only consider 4s epochs without EEG artefacts that are in between 2 other artefact-free expochs of the same sleep-wake state}
          BEGIN
            CASE spec4.state OF
              'w' : state:=1;
              'n' : state:=2;
              'r' : state:=3;
            END;
            inc(c[state]);
            FOR hz:=3 TO 360 DO spec_state[state,hz]:=spec_state[state,hz]+spec4.bin[hz];
          END;
        END;
      END;
      close(EEGf);

      {REFERENCE; these individual references are added to the directory file so that they can be used for other EEG phenotypes without having to recalculate them}
      FOR state:=1 TO 3 DO
      BEGIN
        ref[state]:=0.0;
        FOR hz:=3 to 360 DO
        BEGIN
          spec_state[state,hz]:=spec_state[state,hz]*1000000/c[state];   {average spectra per state in microvolt}
          IF (hz<180) OR (hz>220) THEN ref[state]:=ref[state]+spec_state[state,hz]; {exclude 50Hz noise from reference; 45-55Hz}
        END;
      END;
      reft:=0.0; FOR state:=1 TO 3 DO reft:=reft+(ref[state]*rat[state]); 
	  {determine the weighted average total EEG power in baseline. Note that the same reference is used for all 3 states 
	  so that relative power can still be compared among states}


      {write output for 0.75-45Hz}
      FOR hz:=3 to 180 DO
      BEGIN
        hzr:=hz; hzr:=hzr/4;
        write(OUTf,str:3,m:4,hzr:8:2);
        FOR state:=1 TO 3 DO
        BEGIN
          spec_state[state,hz]:=spec_state[state,hz]*100/reft;
          write(OUTf,spec_state[state,hz]:10:4);
        END;
        writeln(OUTf);
      END;

    END; {mouse}
  END;  {strain}
  close(DIRf);
  close(OUTf);
END.