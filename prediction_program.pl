###########PROGRAM TO PREDICT PROMOTER SEQ##################
print "Enter file name:";
$filename=<STDIN>;
chomp($filename);

open(OUT,">result.txt");

open(IN,"$filename");
@seq=<IN>;
chomp(@seq);

open(RG,"reg_elements_wp_143.txt");
@ele=<RG>;
chomp(@ele);

open(OLFO,"olig_4_result.txt");
@olfour=<OLFO>;
chomp(@olfour);

open(OLFI,"olig_5_result.txt");
@olfive=<OLFI>;
chomp(@olfive);

open(OLS,"olig_6_result.txt");
@olsix=<OLS>;
chomp(@olsix);

###########COUNT###########
$tp=0;
$fn=0;
$tn=0;
$fp=0;

####### AT CONTENT ################

for($i=0;$i<=$#seq;$i++)
{ 
	$AT=0;
	$seq=$seq[$i];
	$nA=$seq=~s/A/A/g;
	$nT=$seq=~s/T/T/g;
	$L=length($seq);
	$AT=($nA+$nT)/$L;
	
	#######  Oligomer 2 #############
	$Two=0;
	%olig_two=qw(AT 1.15105067721062 AG 0.916320380763837 AC 1.07719546516546 AA 1.57382382495911 TA 0.933389704923983 TT -0.845372067456051 TG 0.748898849628016 TC 1.06327883640544 CG 0.9061380106665357 CT 0.910179929065807 CA 1.35610192739895 CC 1.2884985478438 GA 0.924492211188852 GG 0.776020296444568 GC 0.972261649878195 GT 0.678829920032253);
	for($j=0;$j<$L-1;$j++)
	{
		$di=substr($seq,$j,2);
		$Two=$Two+$olig_two{$di};
	}
	#######  Oligomer 3 #############
	$Three=0;
	%olig_three=qw(AAA	2.19357087052682	AAT	1.42938489965253	AAC	1.3243988408165	AAG	1.11654282484211	ATA	1.3358371940336	ATT	1.0814013254552	ATC	1.42455038974512	ATG	0.830692138064037	ACA	1.4348314848954	ACT	0.855583466515688	ACC	1.07126857956764	ACG	0.978440597733758	AGA	1.16339594825659	AGT	0.698220670158809	AGC	1.08631681963254	AGG	0.726516032140645	TAA	1.26859782514878	TAT	1.04536828047266	TAC	0.70837855096536	TAG	0.634625499585839	TTA	0.922622184542767	TTT	0.874793091564541	TTC	0.861213727852337	TTG	0.701775812090807	TCA	1.22207167716054	TCT	0.968270049345717	TCC	1.23724516160148	TCG	0.871332757521293	TGA	0.7836967282558	TGT	0.635915319551454	TGC	0.828992014257806	TGG	0.816808151913257	CAA	1.53931067231775	CAT	1.27378678551217	CAC	1.43048281892299	CAG	1.11081842366074	CTA	0.782300049831612	CTT	0.832039943466391	CTC	1.25778084778938	CTG	0.75772820578232	CCA	1.66964503605642	CCT	0.930962549092423	CCC	1.75248673068352	CCG	1.03345749231896	CGA	0.870666204657128	CGT	0.928671035057955	CGC	1.08105495584439	CGG	0.759964130002622	GAA	1.0918884821573	GAT	0.82166816515665	GAC	0.874393839424546	GAG	0.871979538264791	GTA	0.640743615536427	GTT	0.580219195300044	GTC	0.874836008189901	GTG	0.722746380493772	GCA	1.10247022619325	GCT	0.829737836532412	GCC	1.2512939008731	GCG	0.750315991708378	GGA	0.830904248489047	GGT	0.611598022667714	GGC	0.909024489747885	GGG	0.79786483894624);
	for($j=0;$j<$L-2;$j++)
	{
		$di=substr($seq,$j,3);
		$Three=$Three+$olig_three{$di};
	}
	
	#######  Oligomer 4 #############
	$Four=0;
	for($ol=0;$ol<=$#olfour;$ol++)
	{
		($ola,$olb)=split('\t',$olfour[$ol]);
		$olig_four{$ola}=$olb;
	}
		
	for($j=0;$j<$L-3;$j++)
	{
		$di=substr($seq,$j,4);
		$Four=$Four+$olig_four{$di};
	}
	
	#######  Oligomer 5 #############
	$Five=0;
	for($ol=0;$ol<=$#olfive;$ol++)
	{
		($ola,$olb)=split('\t',$olfive[$ol]);
		$olig_five{$ola}=$olb;
	}
		
	for($j=0;$j<$L-4;$j++)
	{
		$di=substr($seq,$j,5);
		$Five=$Five+$olig_five{$di};
	}
	#######  Oligomer 6 #############
	$Six=0;
	for($ol=0;$ol<=$#olsix;$ol++)
	{
		($ola,$olb)=split('\t',$olsix[$ol]);
		$olig_six{$ola}=$olb;
	}
		
	for($j=0;$j<$L-5;$j++)
	{
		$di=substr($seq,$j,6);
		$Six=$Six+$olig_six{$di};
	}
	
	########## Free energy #########
	$FREEENERGY=0;
	%free=qw(AT -0.88 AG -1.30 AC -1.45 AA -1.00 TA -0.58 TT -1.00 TG -1.44 TC -1.28 CG -2.17 CT -1.28 CA -1.45 CC -1.84 GA -1.30 GG -1.84 GC -2.24 GT -1.44 G 0.98 C 0.98 A 1.03 T 1.03);
	for($j=0;$j<$L-1;$j++)
	{
		$di=substr($seq,$j,2);
		$FREEENERGY=$FREEENERGY+$free{$di};
	}
	$ini=substr($seq,0,1);
	$sym=substr($seq,$L-1,1);
	$FREEENERGY=$FREEENERGY+$free{$ini}+$free{$sym};
	
	
	#################Score ########
	%count=();
	$totcount=0;
	$fr=0;
	$rev=$seq;
	$rev=~tr/ATGC/TACG/;
	@sense=($seq,$rev);
	foreach $sense(@sense)
	{
	for($k=0;$k<=$#ele;$k++)
	{
		($motname,$motseq,$w)=split('\t',$ele[$k]);
		chomp($motseq);
		$motcount=$sense=~s/$motseq/$motseq/g;
		$count{$motname}=$motcount*$w if($motcount!=0);
		$totcount=$totcount+$motcount;
	}
	foreach $x(keys %count)
	{
		$fr=$fr+($count{$x}/$totcount);
	}
	}
		$SCORE=($fr/$L)*100;
	
	##########Skew###########
	$g=$seq=~s/G/G/g;
	$c=$seq=~s/C/C/g;
	$a=$seq=~s/A/A/g;
	$t=$seq=~s/T/T/g;
	$cgskew=($c-$g)/($c+$g);
	$acskew=($a+$c-$g-$t)/($a+$c+$g+$t);
	
	
########### property value ####################
		@val=();
		open(DATA,"diNTproperty_9.txt");
		@data=<DATA>;
		chomp(@data);
		for($p=1;$p<=$#data;$p++)
		{
			@prop=split('\t',$data[$p]);
			@di=qw(0 0 AA AC AG AT CA CC CG CT GA GC GG GT TA TC TG TT);	
			for($pr=2;$pr<=17;$pr++)
			{
				$divalue{$di[$pr]}=$prop[$pr];
			}
			$value=0;
			$l=length($seq);
			for($di=0;$di<=$l-1;$di++)
			{
				$sub=substr($seq,$di,2);
				$value=$value+$divalue{$sub};
			}
			$value=$value/$l;
			push(@val,$value);	
		}
		$Twist=$val[0];
		$MinorGrooveSize=$val[1];
		$Mobilitymajorgrooove=$val[2];
		$Rise_rise=$val[3];
		$Adeninecontent=$val[4];
		$Guaninecontent=$val[5];
		$Thyminecontent=$val[6];
		$Direction=$val[7];
		$Flexibility_shift=$val[8];
		
		
		$sn=$i+1;
			
if($Six <= 274.83512)
{
    if($Six <= 249.87058)
	{
      if($MinorGrooveSize <= 3.203187)
		{
			if($AT <= 0.657371)
			{
				print OUT "$sn\tPROMOTER SEQUENCE\n";
				
				
			}
 			elsif($AT > 0.657371)
			{
				print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";
							
			}
		}
		elsif($MinorGrooveSize > 3.203187)
		{
			if($Six <= 235.231945)
			{
      			if($Flexibility_shift <= 6.353147)
				{
					if($Guaninecontent <= 0.458167)
					{
						print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";
												
					}
					elsif($Guaninecontent > 0.458167)
					{
             			if($Six <= 214.523599)
						{
							print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";
														
						}
             			elsif($Six > 214.523599)
						{
							if($acskew <= -0.155378)
							{
                 				if($SCORE <= 1.477469)
							    {
                  					 if($Three <= 231.103569)
									{
										print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";
																			
									}
									elsif($Three > 231.103569)
									{
										print OUT "$sn\tPROMOTER SEQUENCE\n";
																		
									}
								}
                 				elsif($SCORE > 1.477469)
								{
                  					 if($cgskew <= -0.130435)
									{
                     					if($Twist <= 36.421116)
										{
                       					if($Four <= 224.468504)
											{
                         						if($Four <= 224.123486)
												{
													print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";
																										
												}
                         						elsif($Four > 224.123486)
												{
													print OUT "$sn\tPROMOTER SEQUENCE\n";
													
													
												}
											}
                      						 elsif($Four > 224.468504)
											{
												print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";
																								
											}
										}
										elsif($Twist > 36.421116)
										{
											if($FREEENERGY <= -321.02)
											{
												if($acskew <= -0.163347)
												{
													print OUT "$sn\tPROMOTER SEQUENCE\n";
													
													
												}
												elsif($acskew > -0.163347)
												{
													print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";
																										
												}
											}
											elsif($FREEENERGY > -321.02)
											{
												print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";
																								
											}
										}
									}
									elsif($cgskew > -0.130435)
									{
										if($cgskew <= -0.086207)
										{
											print OUT "$sn\tPROMOTER SEQUENCE\n";
											 
											
										}
										elsif($cgskew > -0.086207)
										{
											print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";
																						
										}
									}
								}
							}
							elsif($acskew > -0.155378)
							{
								if($Flexibility_shift <= 6.071912)
								{
									if($Six <= 233.814536)
									{
										if($Four <= 234.85285)
										{
											if($Rise_rise <= 7.675387)
											{
												if($Six <= 227.976007)
												{
													print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";
														
												}
												elsif($Six > 227.976007)
												{
													print OUT "$sn\tPROMOTER SEQUENCE\n";
													
													
												}
											}
											elsif($Rise_rise > 7.675387)
											{
												print OUT "$sn\tPROMOTER SEQUENCE\n";
												
												
											}
										}
										elsif($Four > 234.85285)
										{
											print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";
																						
										}
									}
									elsif($Six > 233.814536)
									{
										print OUT "$sn\tPROMOTER SEQUENCE\n";
										  
										
									}
								}
								elsif($Flexibility_shift > 6.071912)
								{
									print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";
																		
								}
							}
						}
					}
				}
				elsif($Flexibility_shift > 6.353147)
				{
					print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";
										
				}
			}
			elsif($Six > 235.231945)
			{
				if($Four <= 237.918237)
				{
					if($MinorGrooveSize <= 3.315857)
					{
						print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";
						
						
					}
					elsif($MinorGrooveSize > 3.315857)
					{
						if($Flexibility_shift <= 6.376813)
						{
							if($Guaninecontent <= 0.390438)
							{
								print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";
								 
								
							}
							elsif($Guaninecontent > 0.390438)
							{
								print OUT "$sn\tPROMOTER SEQUENCE\n";
								
								
							}
						}
						elsif($Flexibility_shift > 6.376813)
						{
							if($Three <= 233.815756)
							{
							if($acskew <= -0.179283)
								{
									print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";
																	
								}
							elsif($acskew > -0.179283)
								{
									print OUT "$sn\tPROMOTER SEQUENCE\n";
									
									
								}
							}
							elsif($Three > 233.815756)
							{
							if($SCORE <= 7.473447)
								{
									print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";
																		
								}
								elsif($SCORE > 7.473447)
								{
								if($FREEENERGY <= -386.2)
									{
										print OUT "$sn\tPROMOTER SEQUENCE\n";
										
										
									}
									elsif($FREEENERGY > -386.2)
									{
										print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
										
									}
								}
							}
						}
					}
				}
				elsif($Four > 237.918237)
				{
					if($acskew <= -0.048)
					{
						if($Guaninecontent <= 0.410359)
						{
							if($Thyminecontent <= 0.721116)
							{
								if($Six <= 246.913957)
								{
									if($Four <= 241.601701)
									{
										if($Four <= 241.264028)
										{
											print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";
																															
										}
										elsif($Four > 241.264028)
										{
											print OUT "$sn\tPROMOTER SEQUENCE\n";
											
										}
									}
									elsif($Four > 241.601701)
									{
										print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
										
									}
								}
								elsif($Six > 246.913957)
								{
									if($Twist <= 36.228088)
									{	
										print OUT "$sn\tNOT A PROMOTER SEQUENCE\n"; 
										
									}
									elsif($Twist > 36.228088)
									{
										if($Thyminecontent <= 0.673307)
										{
											print OUT "$sn\tPROMOTER SEQUENCE\n";
																						
										}
										elsif($Thyminecontent > 0.673307)
										{
											if($Six <= 247.231962)
											{
												print OUT "$sn\tPROMOTER SEQUENCE\n";
												
											}
											elsif($Six > 247.231962)
											{
												if($cgskew <= -0.134021)
												{
													if($Five <= 250.002532)
													{
														print OUT "$sn\tPROMOTER SEQUENCE\n";
														
													}
													elsif($Five > 250.002532)
													{
														print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
														
													}
												}
												elsif($cgskew > -0.134021)
												{
													print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";
													
												}
											}
										}
									}
								}
							}
							elsif($Thyminecontent > 0.721116)
							{
								print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";
								
							}
						}
						elsif($Guaninecontent > 0.410359)
						{
							if($Twist <= 36.265697)
							{
								if($acskew <= -0.091633)
								{
									if($Direction <= 17.896414)
									{
										print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";
										
									}
									elsif($Direction > 17.896414)
									{
										if($Guaninecontent <= 0.466135)
										{
											print OUT "$sn\tNOT A PROMOTER SEQUENCE\n"; 
											
										}
										elsif($Guaninecontent > 0.466135)
										{
											if($acskew <= -0.099602)
											{
												print OUT "$sn\tPROMOTER SEQUENCE\n";
												
											}
											elsif($acskew > -0.099602)
											{
												print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";
																														
											}
										}
									}
								}
								elsif($acskew > -0.091633)
								{
									print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";
									
								}
							}
							elsif($Twist > 36.265697)
							{
								if($FREEENERGY <= -311.57)
								{
									if($Flexibility_shift <= 6.229004)
									{
										if($Guaninecontent <= 0.458167)
										{
											print OUT "$sn\tNOT A PROMOTER SEQUENCE\n"; 
											
										}
										elsif($Guaninecontent > 0.458167)
										{
											if($acskew <= -0.195219)
											{
												if($cgskew <= -0.266055)
												{
													print OUT "$sn\tPROMOTER SEQUENCE\n";
													
												}
												elsif($cgskew > -0.266055)
												{
													print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";
													
												}
											}
											elsif($acskew > -0.195219)
											{
												print OUT "$sn\tPROMOTER SEQUENCE\n"; 
												
											}
										}
									}
									elsif($Flexibility_shift > 6.229004)
									{
										if($Six <= 242.10771)
										{
											print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";
											
										}
										elsif($Six > 242.10771)
										{
											if($Adeninecontent <= 0.545817)
											{
												print OUT "$sn\tPROMOTER SEQUENCE\n";
												
											}
											elsif($Adeninecontent > 0.545817)
											{
												print OUT "$sn\tNOT A PROMOTER SEQUENCE\n"; 
												
											}
										}
									}
								}
								elsif($FREEENERGY > -311.57)
								{
									print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";
									
								}
							}
						}
					}
					elsif($acskew > -0.048)
					{
						print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";
						
					}
				}
			}
		}
		
	}
	elsif($Six > 249.87058)
	{
		if($Guaninecontent <= 0.422311)
		{
			if($acskew <= 0.043825)
			{
				if($MinorGrooveSize <= 3.319841)
				{
					if($Flexibility_shift <= 5.800677)
					{
						if($Six <= 263.135039)
						{
							if($Guaninecontent <= 0.394422)
							{
								print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";
								
							}
							elsif($Guaninecontent > 0.394422)
							{
								if($Guaninecontent <= 0.398406)
								{
									print OUT "$sn\tPROMOTER SEQUENCE\n";
									
								}
								elsif($Guaninecontent > 0.398406)
								{
									print OUT "$sn\tNOT A PROMOTER SEQUENCE\n"; 
									
								}
							}
						}
						elsif($Six > 263.135039)
						{
						if($Flexibility_shift <= 5.334223)
							{
								print OUT "$sn\tNOT A PROMOTER SEQUENCE\n"; 
								
							}
							elsif($Flexibility_shift > 5.334223)
							{
							if($Direction <= 12.442231)
								{
									print OUT "$sn\tNOT A PROMOTER SEQUENCE\n"; 
									
								}
								elsif($Direction > 12.442231)
								{
								if($Guaninecontent <= 0.346614)
									{
									if($Rise_rise <= 7.616751)
										{
											if($Guaninecontent <= 0.318725)
											{
												if($Two <= 160.267364)
												{
													if($FREEENERGY <= -284.96)
													{
														print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";
														
													}
													elsif($FREEENERGY > -284.96)
													{
														print OUT "$sn\tPROMOTER SEQUENCE\n";
														
													}
												}
												elsif($Two > 160.267364)
												{
													print OUT "$sn\tPROMOTER SEQUENCE\n";
													
												}
											}
											elsif($Guaninecontent > 0.318725)
											{
												print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";
												
											}
										}
										elsif($Rise_rise > 7.616751)
										{
										if($Three <= 258.692185)
											{
												print OUT "$sn\tNOT A PROMOTER SEQUENCE\n"; 
												
											}
											elsif($Three > 258.692185)
											{
												if($Three <= 261.230756)
												{
													if($Six <= 274.062895)
													{
														if($Guaninecontent <= 0.231076)
														{
															print OUT "$sn\tPROMOTER SEQUENCE\n"; 
															
														}
														elsif($Guaninecontent > 0.231076)
														{
														if($SCORE <= 7.392002)
															{
																if($Mobilitymajorgroove <= 1.070876)
																{
																	print OUT "$sn\tNOT A PROMOTER SEQUENCE\n"; 
																	
																}
																elsif($Mobilitymajorgroove > 1.070876)
																{
																	print OUT "$sn\tPROMOTER SEQUENCE\n";
																																
																}
															}
															elsif($SCORE > 7.392002)
															{
																print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";
																														
															}
														}
													}
													elsif($Six > 274.062895)
													{
														print OUT "$sn\tPROMOTER SEQUENCE\n";
														
													}
												}
												elsif($Three > 261.230756)
												{
													print OUT "$sn\tNOT A PROMOTER SEQUENCE\n"; 
													
												}
											}
										}
									}
									elsif($Guaninecontent > 0.346614)
									{
										if($cgskew <= -0.252033)
										{
											print OUT "$sn\tNOT A PROMOTER SEQUENCE\n"; 
											
										}
										elsif($cgskew > -0.252033)
										{
											if($SCORE <= 6.563419)
											{
												if($Twist <= 36.490359)
												{
													print OUT "$sn\tNOT A PROMOTER SEQUENCE\n"; 
													
												}
												elsif($Twist > 36.490359)
												{
													if($$Twist <= 36.664183)
													{
														print OUT "$sn\tPROMOTER SEQUENCE\n";
																											
													}
													elsif($Twist > 36.664183)
													{
														if($Three <= 250.470153)
														{
															print OUT "$sn\tPROMOTER SEQUENCE\n";
															
														}
														elsif($Three > 250.470153)
														{
															print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
															
														}
													}
												}
											}
											elsif($SCORE > 6.563419)
											{
												print OUT "$sn\tPROMOTER SEQUENCE\n";
												
											}
										}
									}
								}
							}
						}
					
					
					}
					elsif($Flexibility_shift > 5.800677)
					{
						if($AT <= 0.645418)
						{
							if($Rise_rise <= 7.595145)
							{
								if($Mobilitymajorgroove <= 1.062669)
								{
									print OUT "$sn\tPROMOTER SEQUENCE\n";
								}
								elsif($Mobilitymajorgroove > 1.062669)
								{
									print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";
									
								}
							}
							elsif($Rise_rise > 7.595145)
							{
								print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";
							}
						}
						elsif($AT > 0.645418)
						{
							if($acskew <= -0.171315)
							{
								if($Mobilitymajorgroove <= 1.059163)
								{
									print OUT "$sn\tPROMOTER SEQUENCE\n";
								}
								elsif($Mobilitymajorgroove > 1.059163)
								{
									print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";
								}
							}
							elsif($acskew > -0.171315)
							{
								if($Direction <= 10.98008)
								{
									print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";
								}
								elsif($Direction > 10.98008)
								{
									if($acskew <= -0.099602)
									{
										if($acskew <= -0.131474)
										{
											if($Six <= 261.409387)
											{
												if($Mobilitymajorgroove <= 1.064143)
												{
													if($Two <= 168.825836)
													{
														print OUT "$sn\tPROMOTER SEQUENCE\n";
													}
													elsif($Two > 168.825836)
													{
														print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";
													}
												}
												elsif($Mobilitymajorgroove > 1.064143)
												{	
													print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";
												}
											}
											elsif($Six > 261.409387)
											{
												print OUT "$sn\tPROMOTER SEQUENCE\n";
											}
										}
										elsif($acskew > -0.131474)
										{
											if($Guaninecontent <= 0.370518)
											{
												if($SCORE <= 7.871412)
												{
													print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";
												}
												elsif($SCORE > 7.871412)
												{
													if($Mobilitymajorgroove <= 1.061275)
													{
														print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";
													}
													elsif($Mobilitymajorgroove > 1.061275)
													{
														print OUT "$sn\tPROMOTER SEQUENCE\n";
													}
												}
											}
											elsif($Guaninecontent > 0.370518)
											{
												if($Thyminecontent <= 0.701195)
												{
													print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";
												}
												elsif($Thyminecontent > 0.701195)
												{
													print OUT "$sn\tPROMOTER SEQUENCE\n";
												}
											}
										}
									}
									elsif($acskew > -0.099602)
									{
										if($Direction <= 21.035857)
										{
											if($Flexibility_shift <= 5.816494)
											{
												if($AT <= 0.665339)
												{
													print OUT "$sn\tPROMOTER SEQUENCE\n";
													
												}
												elsif($AT > 0.665339)
												{
													print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";
													
												}
											}
											elsif($Flexibility_shift > 5.816494)
											{
												if($AT <= 0.653386)
												{												
													print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";
													
												}
											    elsif($AT > 0.653386)
												{
													print OUT "$sn\tPROMOTER SEQUENCE\n";
													
												}
											}
										}
										elsif($Direction > 21.035857)
										{
											print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";
										}
									}	
							    }
							}
						}
					}
					
				}
				elsif($MinorGrooveSize > 3.319841)
				{
					if($Six <= 261.74334)
					{
						if($Direction <= 6.422311)
						{
							if($acskew <= -0.091633)
							{
								if($Rise_rise <= 7.64913)
								{;
									if($cgskew <= 0.02994)
									{
										print OUT "$sn\tPROMOTER SEQUENCE\n";
										
									}
									elsif($cgskew > 0.02994)
									{
										print OUT "$sn\tNOT A PROMOTER SEQUENCE\n"; 
									}
								}
								elsif($Rise_rise > 7.64913)
								{
										print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";
								}
							}
							elsif($acskew > -0.091633)
							{
								print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";
							}
						}
						elsif($Direction > 6.422311)
						{
							if($acskew <= -0.035857)
							{
								if($Flexibility_shift <= 5.793227)
								{
									if($Twist <= 36.739562)
									{
										print OUT "$sn\tNOT A PROMOTER SEQUENCE\n"; 
									}
									elsif($Twist > 36.739562)
									{
										if($Two <= 189.600815)
										{
												print OUT "$sn\tPROMOTER SEQUENCE\n";
												 
										}
										elsif($Two > 189.600815)
										{
												print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";
										}
									}
								}
								elsif($Flexibility_shift > 5.793227)
								{
									if($Thyminecontent <= 0.800797)
									{
										if($acskew <= -0.139442)
										{
											print OUT "$sn\tPROMOTER SEQUENCE\n";
										}
										elsif($acskew > -0.139442)
										{
											if($AT <= 0.63745)
											{
												print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
											}
											elsif($AT > 0.63745)
											{
												print OUT "$sn\tPROMOTER SEQUENCE\n";
											}
										}
									}
									elsif($Thyminecontent > 0.800797)
									{
										if($Rise_rise <= 7.648352)
										{
											print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";
										}
										elsif($Rise_rise > 7.648352)
										{
											print OUT "$sn\tPROMOTER SEQUENCE\n";
										}
									}
								}
							}
							elsif($acskew > -0.035857)
							{
								if($Flexibility_shift <= 6.095179)
								{
									if($FREEENERGY <= -327.94)
									{
										if($acskew <= 0.020576)
										{
											print OUT "$sn\tPROMOTER SEQUENCE\n";
										}
										elsif($acskew > 0.020576)
										{
											print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";
										}
									}
									elsif($FREEENERGY > -327.94)
									{
										print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
									}
								}
								elsif($Flexibility_shift > 6.095179)
								{
										print OUT "$sn\tNOT A PROMOTER SEQUENCE\n"; 
								}
							}
						}
					}
					elsif($Six > 261.74334)
					{
						if($Direction <= 4.031873)
                        {
							print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
						}
						elsif($Direction > 4.031873)
						{
						if($Flexibility_shift <= 5.698127)
							{
								if($AT <= 0.593625)
								{
									print OUT "$sn\tPROMOTER SEQUENCE\n";
								}
								elsif($AT > 0.593625)
								{
									if($AT <= 0.641434)
									{
										if($acskew <= -0.056)
										{
											if($Direction <= 9.988048)
											{
												print OUT "$sn\tNOT A PROMOTER SEQUENCE\n"; 
											}
											elsif($Direction > 9.988048)
											{
												if($SCORE <= 4.209021)
												{
													print OUT "$sn\tPROMOTER SEQUENCE\n";
												}
												elsif($SCORE > 4.209021)
												{
													if($Twist <= 36.695817)
													{
														if($FREEENERGY <= -311.22)
														{
															print OUT "$sn\tNOT A PROMOTER SEQUENCE\n"; 
														}
														elsif($FREEENERGY > -311.22)
														{
															print OUT "$sn\tPROMOTER SEQUENCE\n";
														}
													}
													elsif($Twist > 36.695817)
													{
														print OUT "$sn\tPROMOTER SEQUENCE\n";
													}
												}
											}
											
										}
										elsif($acskew > -0.056)
										{
											print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";
										}
									}
									elsif($AT > 0.641434)
									{
										print OUT "$sn\tPROMOTER SEQUENCE\n";
									}
								}
							}
							elsif($Flexibility_shift > 5.698127)
							{
								if($acskew <= 0.020576)
								{
									if($Six <= 270.806158)
									{
											print OUT "$sn\tPROMOTER SEQUENCE\n"; 
									}
									elsif($Six > 270.806158)
									{
										if($Guaninecontent <= 0.40239)
										{
											if($Direction <= 10.848606)
											{
												if($Five <= 261.619944)
												{
														print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";
												}
												elsif($Five > 261.619944)
												{
													print OUT "$sn\tPROMOTER SEQUENCE\n";
												}
											}
											if($Direction > 10.848606)
											{
												print OUT "$sn\tPROMOTER SEQUENCE\n";  
											}
										}
										elsif($Guaninecontent > 0.40239)
										{
											if($Guaninecontent <= 0.418327)
											{
												print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
											}
											elsif($Guaninecontent > 0.418327)
											{
													print OUT "$sn\tPROMOTER SEQUENCE\n"; 
											}
										}
									}
								}
								elsif($acskew > 0.020576)
								{
									if($Flexibility_shift <= 6.017849)
									{
										if($Thyminecontent <= 0.541833)
										{
											print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
										}
										elsif($Thyminecontent > 0.541833)
										{
											print OUT "$sn\tPROMOTER SEQUENCE\n";  
										}
									}
									elsif($Flexibility_shift > 6.017849)
									{
										print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
									}
								}
							}
						   
						}
					}
				}
			}
			elsif($acskew > 0.043825)
			{
				if($Rise_rise <= 7.632794)
				{
					if($MinorGrooveSize <= 3.331793)
					{
						if($Flexibility_shift <= 5.902669)
						{
							print OUT "$sn\tPROMOTER SEQUENCE\n";
						}
						elsif($Flexibility_shift > 5.902669)
						{
							print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
						}
					}
					elsif($MinorGrooveSize > 3.331793)
					{
						print OUT "$sn\tNOT A PROMOTER SEQUENCE\n"; 
					}
				}
				elsif($Rise_rise > 7.632794)
				{
					if($Rise_rise <= 7.663522)
					{
						if($Five <= 267.226346)
						{
							print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
						}
						elsif($Five > 267.226346)
						{
							if($Five <= 270.915807)
							{
								print OUT "$sn\tPROMOTER SEQUENCE\n";
							}
							elsif($Five > 270.915807)
							{
								print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
							}
						}
					}
					elsif($Rise_rise > 7.663522)
					{
						print OUT "$sn\tPROMOTER SEQUENCE\n"; 
					}
				}
			}
		}
		elsif($Guaninecontent > 0.422311)
		{
			if($acskew <= 0.011952)
			{
				if($Six <= 262.954731)
				{
					if($Three <= 243.725785)
					{
						if($Two <= 230.443862)
						{
							if($FREEENERGY <= -352.96)
							{
								print OUT "$sn\tPROMOTER SEQUENCE\n";  
							}
							elsif($FREEENERGY > -352.96)
							{
								if($Flexibility_shift <= 6.445777)
								{
									if($Two <= 197.853057)
									{
										if($acskew <= -0.155378)
										{
											if($Adeninecontent <= 0.338645)
											{
												print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
											}
											elsif($Adeninecontent > 0.338645)
											{
												print OUT "$sn\tPROMOTER SEQUENCE\n";
											}
										}
										elsif($acskew > -0.155378)
										{
											print OUT "$sn\tNOT A PROMOTER SEQUENCE\n"; 
										}
									}
									elsif($Two > 197.853057)
									{
										print OUT "$sn\tPROMOTER SEQUENCE\n"; 
									}
								}
								elsif($Flexibility_shift > 6.445777)
								{
									print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
								}
							}
						}
						elsif($Two > 230.443862)
						{
							if($Four <= 245.203018)
							{
								print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
							}
							elsif($Four > 245.203018)
							{
								print OUT "$sn\tPROMOTER SEQUENCE\n";
							}
						}
					}
					elsif($Three > 243.725785)
					{
						if($MinorGrooveSize <= 3.590757)
						{
							if($Direction <= 9.541833)
							{
								if($Six <= 255.580506)
								{
									if($acskew <= -0.155378)
									{
										if($AT <= 0.649402)
										{
											print OUT "$sn\tPROMOTER SEQUENCE\n"; 
										}
										elsif($AT > 0.649402)
										{
											print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
										}
									}
									elsif($acskew > -0.155378)
									{
										print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
									}
								}
								elsif($Six > 255.580506)
								{
									if($AT <= 0.450199)
									{
										print OUT "$sn\tPROMOTER SEQUENCE\n";
									}
									elsif($AT > 0.450199)
									{
										if($acskew <= -0.016)
										{
											if($Three <= 255.263639)
											{
												if($Guaninecontent <= 0.49004)
												{
													if($AT <= 0.625498)
													{
														if($Rise_rise <= 7.593627)
														{
															if($Five <= 256.307222)
															{
																if($Four <= 253.08674)
																{
																	print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
																}
																elsif($Four > 253.08674)
																{
																	print OUT "$sn\tPROMOTER SEQUENCE\n";
																}
															}
															elsif($Five > 256.307222)
															{
																print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";
															}
														}
														elsif($Rise_rise > 7.593627)
														{
															print OUT "$sn\tPROMOTER SEQUENCE\n"; 
														}
													}
													elsif($AT > 0.625498)
													{
														print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
													}
												}
												elsif($Guaninecontent > 0.49004)
												{
													print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
												}
											}
											elsif($Three > 255.263639)
											{
												print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
											}
										}
										elsif($acskew > -0.016)
										{
											print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
										}
									}
								}
							}
							elsif($Direction > 9.541833)
							{
								if($Rise_rise <= 7.576664)
								{
									if($Twist <= 35.735458)
									{
										print OUT "$sn\tPROMOTER SEQUENCE\n";
									}
									elsif($Twist > 35.735458)
									{
										if($AT <= 0.386454)
										{
											print OUT "$sn\tPROMOTER SEQUENCE\n";
										}
										elsif($AT > 0.386454)
										{
											if($Direction <= 16.306773)
											{
												if($Six <= 260.442232)
												{
													print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
												}
												elsif($Six > 260.442232)
												{
													print OUT "$sn\tPROMOTER SEQUENCE\n";
												}
											}
											elsif($Direction > 16.306773)
											{
												if($Rise_rise <= 7.513908)
												{
													print OUT "$sn\tNOT A PROMOTER SEQUENCE\n"; 
												}
												elsif($Rise_rise > 7.513908)
												{
												    if($Rise_rise <= 7.565008)
													{
														if($Direction <= 17.756972)
														{
															if($Four <= 249.421985)
															{
																print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
															}
															elsif($Four > 249.421985)
															{
																print OUT "$sn\tPROMOTER SEQUENCE\n";
															}
														}
														elsif($Direction > 17.756972)
														{
															print OUT "$sn\tPROMOTER SEQUENCE\n";  
														}
													}
													elsif($Rise_rise > 7.565008)
													{
														print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
													}
												}
											}
										}
									}
								}
								elsif($Rise_rise > 7.576664)
								{
								if($Flexibility_shift <= 5.535378)
									{
										print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
									}
									elsif($Flexibility_shift > 5.535378)
									{
										if($Two <= 175.371483)
										{
											print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
										}
										elsif($Two > 175.371483)
										{
											if($Three <= 245.764137)
											{
												print OUT "$sn\tPROMOTER SEQUENCE\n"; 
											}
											elsif($Three > 245.764137)
											{
												if($Six <= 255.19317)
												{
													if($Flexibility_shift <= 5.61745)
													{
														print OUT "$sn\tPROMOTER SEQUENCE\n";
													}
													elsif($Flexibility_shift > 5.61745)
													{
														if($Adeninecontent <= 0.486056)
														{
															if($Three <= 247.202184)
															{
																print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
															}
															elsif($Three > 247.202184)
															{
																print OUT "$sn\tPROMOTER SEQUENCE\n"; 
															}
														}
														elsif($Adeninecontent > 0.486056)
														{
															if($Three <= 247.29203)
															{
																print OUT "$sn\tPROMOTER SEQUENCE\n";
															}
															elsif($Three > 247.29203)
															{
																if($Thyminecontent <= 0.561753)
																{
																	print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
																}
																elsif($Thyminecontent > 0.561753)
																{
																	if($Two <= 197.984636)
																	{
																		print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
																	}
																	elsif($Two > 197.984636)
																	{
																		if($Thyminecontent <= 0.641434)
																		{
																			if($Two <= 206.068624)
																			{
																				print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
																			}
																			elsif($Two > 206.068624)
																			{
																				print OUT "$sn\tPROMOTER SEQUENCE\n";
																			}
																		}
																		elsif($Thyminecontent > 0.641434)
																		{
																			print OUT "$sn\tPROMOTER SEQUENCE\n";  
																		}
																	}
																}
															}
														}
													}
												}
												elsif($Six > 255.19317)
												{
													if($FREEENERGY <= -336.78)
													{
														print OUT "$sn\tPROMOTER SEQUENCE\n";
													}
													elsif($FREEENERGY > -336.78)
													{
														if($Mobilitymajorgroove <= 1.053785)
														{
															print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
														}
														elsif($Mobilitymajorgroove > 1.053785)
														{
														if($Direction <= 14.717131)
															{
																if($Direction <= 12.358566)
																{
																	print OUT "$sn\tPROMOTER SEQUENCE\n"; 
																}
																elsif($Direction > 12.358566)
																{
																	if($Guaninecontent <= 0.434263)
																	{
																		print OUT "$sn\tPROMOTER SEQUENCE\n"; 
																	}
																	elsif($Guaninecontent > 0.434263)
																	{
																		print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
																	}
																}
															}
															elsif($Direction > 14.717131)
															{
																print OUT "$sn\tPROMOTER SEQUENCE\n"; 
															}
														}
													}
												}
											}
										}
									}
								}
							}
						}
						elsif($MinorGrooveSize > 3.590757)
						{
							print OUT "$sn\tNOT A PROMOTER SEQUENCE\n"; 
						}
					}
				}
				elsif($Six > 262.954731)
				{
					if($Direction <= 8.912351)
					{
						if($Direction <= -4.816733)
						{
							print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
						}
						elsif($Direction > -4.816733)
						{
							if($Flexibility_shift <= 6.168327)
							{
								if($Guaninecontent <= 0.426295)
								{
									print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
								}
								elsif($Guaninecontent > 0.426295)
								{
									print OUT "$sn\tPROMOTER SEQUENCE\n";
								}
							}
							elsif($Flexibility_shift > 6.168327)
							{
								if($Adeninecontent <= 0.426295)
								{
									if($Mobilitymajorgroove <= 1.033187)
									{
										print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
									}
									elsif($Mobilitymajorgroove > 1.033187)
									{
										print OUT "$sn\tPROMOTER SEQUENCE\n";
									}
								}
								elsif($Adeninecontent > 0.426295)
								{
									if($Twist <= 35.9051)
									{
										if($AT <= 0.537849)
										{
											print OUT "$sn\tPROMOTER SEQUENCE\n";
										}
										elsif($AT > 0.537849)
										{
											print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";
										}
									}
									elsif($Twist > 35.9051)
									{
										print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";
									}
								}
							}
						}
					}
					elsif($Direction > 8.912351)
					{
						if($Flexibility_shift <= 6.410598)
						{
							if($FREEENERGY <= -314.26)
							{
								print OUT "$sn\tPROMOTER SEQUENCE\n"; 
							}
							elsif($FREEENERGY > -314.26)
							{
								if($acskew <= -0.155378)
								{
									print OUT "$sn\tPROMOTER SEQUENCE\n";
								}
								elsif($acskew > -0.155378)
								{
									if($Three <= 256.318514)
									{
										print OUT "$sn\tNOT A PROMOTER SEQUENCE\n"; 
									}
										elsif($Three > 256.318514)
									{
										if($Flexibility_shift <= 5.750558)
										{
											print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";
										}
										elsif($Flexibility_shift > 5.750558)
										{
											print OUT "$sn\tPROMOTER SEQUENCE\n";
										}
									}
								}
							}
						}
						elsif($Flexibility_shift > 6.410598)
						{
							if($Twist <= 35.859681)
							{
								print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";
							}
							elsif($Twist > 35.859681)
							{
								if($cgskew <= -0.095238)
								{
									print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";
								}
								elsif($cgskew > -0.095238)
								{
									if($acskew <= -0.011952)
									{
										print OUT "$sn\tPROMOTER SEQUENCE\n";
									}
									elsif($acskew > -0.011952)
									{
										print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
									}
								}
							}
						}
					}
				}
			}
			elsif($acskew > 0.011952)
			{
				if($Rise_rise <= 7.609185)
				{
					if($Six <= 260.316799)
					{
						print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
					}
					elsif($Six > 260.316799)
					{
						if($FREEENERGY <= -350.56)
						{
							if($Flexibility_shift <= 6.47506)
							{
								if($Four <= 261.998736)
								{
									if($Flexibility_shift <= 6.326255)
									{
										print OUT "$sn\tPROMOTER SEQUENCE\n";
									}
									elsif($Flexibility_shift > 6.326255)
									{
										if($Guaninecontent <= 0.466135)
										{
											print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";
										}
										elsif($Guaninecontent > 0.466135)
										{
											if($AT <= 0.278884)
											{
												print OUT "$sn\tNOT A PROMOTER SEQUENCE\n"; 
											}
											elsif($AT > 0.278884)
											{
												if($Three <= 258.72332)
												{
													if($Direction <= 12.856574)
													{
														if($Five <= 259.948287)
														{
															print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";
														}
														if($Five > 259.948287)
														{
															if($Mobilitymajorgroove <= 1.042988)
															{
																print OUT "$sn\tPROMOTER SEQUENCE\n";
															}
															elsif($Mobilitymajorgroove > 1.042988)
															{
																print OUT "$sn\tNOT A PROMOTER SEQUENCE\n"; 
															}
														}
													}
													elsif($Direction > 12.856574)
													{
														print OUT "$sn\tPROMOTER SEQUENCE\n";
													}
												}
												elsif($Three > 258.72332)
												{
													print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
												}
											}
										}
									}
								}
								elsif($Four > 261.998736)
								{
									print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";;  
								}
							}
							elsif($Flexibility_shift > 6.47506)
							{
								if($Mobilitymajorgroove <= 1.02753)
								{
									print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
								}
								elsif($Mobilitymajorgroove > 1.02753)
								{
								if($SCORE <= 1.836536)
									{
									if($Two <= 242.411008)
										{
											if($acskew <= 0.059761)
											{
												print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
											}
											elsif($acskew > 0.059761)
											{
												print OUT "$sn\tPROMOTER SEQUENCE\n";
											}
										}
										elsif($Two > 242.411008)
										{
											print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
										}
									}
									elsif($SCORE > 1.836536)
									{
										print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
									}
								}
							}
						}
						elsif($FREEENERGY > -350.56)
						{
							if($Rise_rise <= 7.581287)
							{
								print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
							}
							elsif($Rise_rise > 7.581287)
							{
								if($Flexibility_shift <= 6.240438)
								{
									if($acskew <= 0.043825)
									{
										print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
									}
									elsif($acskew > 0.043825)
								    {
									 print OUT "$sn\tPROMOTER SEQUENCE\n";
								    }
								}
								elsif($Flexibility_shift > 6.240438)
								{
									print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
								}
							}
							
						}
					}
				}
				elsif($Rise_rise > 7.609185)
				{
					if($Flexibility_shift <= 6.006135)
					{
						print OUT "$sn\tPROMOTER SEQUENCE\n";
					}
					elsif($Flexibility_shift > 6.006135)
					{
						if($Adeninecontent <= 0.621514)
						{
							if($Six <= 266.2604)
							{
								if($Mobilitymajorgroove <= 1.052231)
								{
									if($Twist <= 36.151474)
									{
										print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
									}
									elsif($Twist > 36.151474)
									{
										print OUT "$sn\tPROMOTER SEQUENCE\n";
									}
								}
								elsif($Mobilitymajorgroove > 1.052231)
								{
									print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
								}
							}
							elsif($Six > 266.2604)
							{
								print OUT "$sn\tPROMOTER SEQUENCE\n";  
							}
						}
						elsif($Adeninecontent > 0.621514)
						{
							print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
						}
					}
				}
			}
		}
	}	
}
elsif($Six>274.83512)
{
	if($Six <= 302.744367)
	{
		if($MinorGrooveSize <= 3.29992)
		{
			if($MinorGrooveSize <= 3.252112)
			{	
				if($AT <= 0.752988)
				{
					if($Rise_rise <= 7.62465)
					{
						if($cgskew <= 0.233083)
						{
							if($acskew <= -0.067729)
							{
								if($FREEENERGY <= -285.2)
								{
									if($Direction <= 11.976096)
									{
										print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
									}
									elsif($Direction > 11.976096)
									{
										print OUT "$sn\tPROMOTER SEQUENCE\n"; 
									}
								}
								elsif($FREEENERGY > -285.2)
								{
									print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
								}
							}
							elsif($acskew > -0.067729)
							{
								print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
							}
						}
						elsif($cgskew > 0.233083)
						{
							print OUT "$sn\tPROMOTER SEQUENCE\n";  
						}
					}
					elsif($Rise_rise > 7.62465)
					{
						print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
					}
				}
				elsif($AT > 0.752988)
				{
					print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
				}
			}
			elsif($MinorGrooveSize > 3.252112)
			{
				if($Rise_rise <= 7.696737)
				{
					if($Six <= 286.421135)
					{
						if($MinorGrooveSize <= 3.28)
						{
							if($AT <= 0.685259)
							{
								print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
							}
							elsif($AT > 0.685259)
							{
								if($cgskew <= -0.013699)
								{
									if($Rise_rise <= 7.662964)
									{
										if($Direction <= 6.685259)
										{
											print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
										}
										elsif($Direction > 6.685259)
										{
										if($cgskew <= -0.027523)
											{
												if($Five <= 270.711008)
												{
													if($SCORE <= 7.90052)
													{
														print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
													}
													elsif($SCORE > 7.90052)
													{
														print OUT "$sn\tPROMOTER SEQUENCE\n"; 
													}
												}
												elsif($Five > 270.711008)
												{
													print OUT "$sn\tPROMOTER SEQUENCE\n";
												}
											}
											elsif($cgskew > -0.027523)
											{
												print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
											}
										}
									}
									elsif($Rise_rise > 7.662964)
									{
										if($acskew <= -0.016)
										{
											print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
										}
										elsif($acskew > -0.016)
										{
											print OUT "$sn\tPROMOTER SEQUENCE\n";  
										}
									}
								}
								elsif($cgskew > -0.013699)
								{
									if($MinorGrooveSize <= 3.276016)
									{
										if($Four <= 265.008046)
										{
											print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
										}
										elsif($Four > 265.008046)
										{
											if($Four <= 271.749302)
											{
												if($Flexibility_shift <= 5.412629)
												{
													print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
												}
												elsif($Flexibility_shift > 5.412629)
												{
													if($AT <= 0.697211)
													{
														if($Five <= 273.063652)
														{
															print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
														}
														elsif($Five > 273.063652)
														{
															if($MinorGrooveSize <= 3.272032)
															{
																if($AT <= 0.693227)
																{
																	if($Three <= 260.827802)
																	{
																		print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
																	}
																	elsif($ Three > 260.827802)
																	{
																		print OUT "$sn\tPROMOTER SEQUENCE\n";  
																	}
																}
																elsif($AT > 0.693227)
																{
																	if($Three <= 264.422539)
																	{
																		print OUT "$sn\tPROMOTER SEQUENCE\n";
																	}
																	elsif($Three > 264.422539)
																	{
																		print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
																	}
																}
															}
															elsif($MinorGrooveSize > 3.272032)
															{
																print OUT "$sn\tPROMOTER SEQUENCE\n";
															}
														}
													}
													elsif($AT > 0.697211)
													{
														print OUT "$sn\tPROMOTER SEQUENCE\n";
													}
												}
											}
											elsif($Four > 271.749302)
											{
												if($AT <= 0.709163)
												{
													if($MinorGrooveSize <= 3.272032)
													{
														print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
													}
													elsif($MinorGrooveSize > 3.272032)
													{
														print OUT "$sn\tPROMOTER SEQUENCE\n"; 
													}
												}
												elsif($AT > 0.709163)
												{
													print OUT "$sn\tPROMOTER SEQUENCE\n"; 
												}
											}
										}
									}
									elsif($MinorGrooveSize > 3.276016)
									{
										print OUT "$sn\tPROMOTER SEQUENCE\n";  
									}
								}
							}
						}
						elsif($MinorGrooveSize > 3.28)
						{
							if($AT <= 0.681275)
							{
								if($Four <= 270.429121)
								{
									if($Six <= 279.197047)	
									{
										if($Three <= 254.915926)
										{
											print OUT "$sn\tPROMOTER SEQUENCE\n";
											 
										}
										elsif($Three > 254.915926)
										{
											if($Thyminecontent <= 0.756972)
											{
												if($Adeninecontent <= 0.613546)
												{
													print OUT "$sn\tPROMOTER SEQUENCE\n"; 
												}
												elsif($Adeninecontent > 0.613546)
												{
													if($cgskew <= -0.042424)
													{
														if($AT <= 0.665339)
														{
															print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
														}
														elsif($AT > 0.665339)
														{
															print OUT "$sn\tPROMOTER SEQUENCE\n"; 
														}
													}
													elsif($cgskew > -0.042424)
													{
														if($AT <= 0.677291)
														{
															print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
														}
														elsif($AT > 0.677291)
														{
															print OUT "$sn\tPROMOTER SEQUENCE\n"; 
														}
													}
												}
											}
											elsif($Thyminecontent > 0.756972)
											{
												print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
											}
										}
									}
									elsif($Six > 279.197047)
									{
										if($Guaninecontent <= 0.278884)
										{
											if($Direction <= 7.613546)
											{
												if($cgskew <= 0.3)
												{
													print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
												}
												elsif($cgskew > 0.3)
												{
													print OUT "$sn\tPROMOTER SEQUENCE\n";  
												}
											}
											elsif($Direction > 7.613546)
											{
												print OUT "$sn\tPROMOTER SEQUENCE\n";   
											}
										}
										elsif($Guaninecontent > 0.278884)
										{
											if($Two <= 191.539535)
											{
												if($MinorGrooveSize <= 3.295936)
												{
													print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
												}
												elsif($MinorGrooveSize > 3.295936)
												{
													if($Adeninecontent <= 0.61753)
													{
														print OUT "$sn\tPROMOTER SEQUENCE\n";   
													}
													elsif($Adeninecontent > 0.61753)
													{
														print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
													}
												}
											}
											elsif($Two > 191.539535)
											{
												print OUT "$sn\tPROMOTER SEQUENCE\n";   
											}
											
										}
									}
								}
								elsif($Four > 270.429121)
								{
									if($Six <= 282.541198)
									{
										print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
									}
									elsif($Six > 282.541198)
									{
										if($FREEENERGY <= -304.4)
										{
											print OUT "$sn\tPROMOTER SEQUENCE\n";   
										}
										elsif($FREEENERGY > -304.4)
										{
											if($Twist <= 36.497769)
											{
												print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
											}
											elsif($Twist > 36.497769)
											{
												if($Direction <= 6.055777)
												{
													print OUT "$sn\tPROMOTER SEQUENCE\n";  
												}
												elsif($Direction > 6.055777)
												{
													print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
												}
											}
										}
									}
								}
							}
							elsif($AT > 0.681275)
							{
								print OUT "$sn\tPROMOTER SEQUENCE\n";  
							}
						}
					}
					elsif($Six > 286.421135)
					{
						if($AT <= 0.665339)
						{
							if($Rise_rise <= 7.585959)
							{
								print OUT "$sn\tPROMOTER SEQUENCE\n"; 
							}
							elsif($Rise_rise > 7.585959)
							{
								if($Guaninecontent <= 0.282869)
								{
								if($acskew <= 0.020576)
									{
										print OUT "$sn\tPROMOTER SEQUENCE\n";   
									}
									elsif($acskew > 0.020576)
									{
										print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
									}
								}
								elsif($Guaninecontent > 0.282869)
								{
									print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
								}
							}
						}
						elsif($AT > 0.665339)
						{
							if($MinorGrooveSize <= 3.290837)
							{
								if($AT <= 0.697211)
								{
									if($MinorGrooveSize <= 3.268048)
									{
										if($Mobilitymajorgroove <= 1.070438)
										{
											if($Flexibility_shift <= 5.6749)
											{
												if($acskew <= -0.011952)
												{
													print OUT "$sn\tNOT A PROMOTER SEQUENCE\n"; 
												}
												elsif($acskew > -0.011952)
												{
													print OUT "$sn\tPROMOTER SEQUENCE\n"; 
												}
											}
											elsif($Flexibility_shift > 5.6749)
											{
												print OUT "$sn\tPROMOTER SEQUENCE\n";   
											}
										}
										elsif($Mobilitymajorgroove > 1.070438)
										{
											print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
										}
									}
									elsif($MinorGrooveSize > 3.268048)
									{
										if($AT <= 0.693227)
										{
											if($MinorGrooveSize <= 3.270916)
											{
												if($FREEENERGY <= -296.9)
												{
													print OUT "$sn\tPROMOTER SEQUENCE\n";  
												}
												elsif($FREEENERGY > -296.9)
												{
													print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
												}
											}
											elsif($MinorGrooveSize > 3.270916)
											{
												if($Adeninecontent <= 0.513944)
												{
													if($Twist <= 36.619203)
													{
														print OUT "$sn\tPROMOTER SEQUENCE\n"; 
													}
													elsif($Twist > 36.619203)
													{
														print OUT "$sn\tNOT A PROMOTER SEQUENCE\n"; 
													}
												}
												elsif($Adeninecontent > 0.513944)
												{
													if($cgskew <= 0.094017)
													{
														if($Five <= 278.364073)
														{
															print OUT "$sn\tPROMOTER SEQUENCE\n";  
														}
														elsif($Five > 278.364073)
														{
															if($AT <= 0.685259)
															{
																if($Guaninecontent <= 0.306773)
																{
																	print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
																}
																elsif($Guaninecontent > 0.306773)
																{
																	if($Rise_rise <= 7.683025)
																	{
																		if($AT <= 0.673307)
																		{
																			if($Flexibility_shift <= 5.667052)
																			{
																				print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
																			}
																			elsif($Flexibility_shift > 5.667052)
																			{
																				if($FREEENERGY <= -301.6)
																				{
																					print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
																				}
																				elsif($FREEENERGY > -301.6)
																				{
																					print OUT "$sn\tPROMOTER SEQUENCE\n";   
																				}
																			}
																		}
																		elsif($AT > 0.673307)
																		{
																			if($MinorGrooveSize <= 3.28)
																			{
																				print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
																			}
																			elsif($MinorGrooveSize > 3.28)
																			{
																				if($Direction <= 5.258964)
																				{
																					if($Six <= 289.146556)
																					{
																						print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
																					}
																					elsif($Six > 289.146556)
																					{
																						if($Four <= 278.551671)
																						{
																							print OUT "$sn\tPROMOTER SEQUENCE\n";   
																						}
																						elsif($Four > 278.551671)
																						{
																							print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
																						}
																					}
																				}
																				elsif($Direction > 5.258964)
																				{
																					print OUT "$sn\tPROMOTER SEQUENCE\n";   
																				}
																			}
																		}
																	}
																	elsif($Rise_rise > 7.683025)
																	{
																		print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
																	}
																}
															}
															elsif($AT > 0.685259)
															{
																if($MinorGrooveSize <= 3.276016)
																{
																	if($AT <= 0.689243)
																	{
																		if($Five <= 285.308983)
																		{
																			print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
																		}
																		elsif($Five > 285.308983)
																		{
																			print OUT "$sn\tPROMOTER SEQUENCE\n";  
																		}
																	}
																	elsif($AT > 0.689243)
																	{
																		if($Guaninecontent <= 0.2749)
																		{
																			print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
																		}
																		elsif($Guaninecontent > 0.2749)
																		{
																			print OUT "$sn\tPROMOTER SEQUENCE\n";   
																		}
																	}
																}
																elsif($MinorGrooveSize > 3.276016)
																{
																	print OUT "$sn\tPROMOTER SEQUENCE\n"; 
																}
															}
														}
													}
													elsif($cgskew > 0.094017)
													{
														if($Thyminecontent <= 0.713147)
														{
															if($Six <= 300.580321)
															{
																if($Five <= 289.138469)
																{
																	if($Two <= 211.512193)
																	{
																		print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
																	}
																	elsif($Two > 211.512193)
																	{
																		if($Flexibility_shift <= 5.587689)
																		{
																			print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
																		}
																		elsif($Flexibility_shift > 5.587689)
																		{
																			print OUT "$sn\tPROMOTER SEQUENCE\n";  
																		}
																	}
																}
																elsif($Five > 289.138469)
																{
																	print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
																}
															}
															elsif($Six > 300.580321)
															{
																print OUT "$sn\tPROMOTER SEQUENCE\n";  
															}
														}
														elsif($Thyminecontent > 0.713147)
														{
															if($Two <= 184.165465)
															{
																if($Guaninecontent <= 0.266932)
																{
																	print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
																}
																elsif($Guaninecontent > 0.266932)
																{
																	print OUT "$sn\tPROMOTER SEQUENCE\n";  
																}
															}
															elsif($Two > 184.165465)
															{
																print OUT "$sn\tPROMOTER SEQUENCE\n";  
															}
														}
													}
												}
											}
										}
										elsif($AT > 0.693227)
										{
											print OUT "$sn\tPROMOTER SEQUENCE\n";  
										}
									}
								}
								elsif($AT > 0.697211)
								{
									if($MinorGrooveSize <= 3.264064)
									{
										if($AT <= 0.709163)
										{
											if($MinorGrooveSize <= 3.256096)
											{
												if($SCORE <= 8.054256)
												{
													print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
												}
												elsif($SCORE > 8.054256)
												{
													if($Six <= 297.716261)
													{
														print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
													}
													elsif($Six > 297.716261)
													{
														print OUT "$sn\tPROMOTER SEQUENCE\n";  
													}
												}
											}
											elsif($MinorGrooveSize > 3.256096)
											{
												if($AT <= 0.705179)
												{
													if($MinorGrooveSize <= 3.26008)
													{
														if($cgskew <= 0.135135)
														{
															if($Direction <= 8.099602)
															{
																print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
															}
															elsif($Direction > 8.099602)
															{
																print OUT "$sn\tPROMOTER SEQUENCE\n"; 
															}
														}
														elsif($cgskew > 0.135135)
														{
															print OUT "$sn\tPROMOTER SEQUENCE\n";   
														}
													}
													elsif($MinorGrooveSize > 3.26008)
													{
														if($AT <= 0.701195)
														{
															if($SCORE <= 7.145959)
															{
																print OUT "$sn\tPROMOTER SEQUENCE\n";   
															}
															elsif($SCORE > 7.145959)
															{
																print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
															}
														}
														elsif($AT > 0.701195)
														{
															print OUT "$sn\tPROMOTER SEQUENCE\n";   
														}
													}
												}
												elsif($AT > 0.705179)
												{
													print OUT "$sn\tPROMOTER SEQUENCE\n";   
												}
											}
										}
										elsif($AT > 0.709163)
										{
											print OUT "$sn\tPROMOTER SEQUENCE\n";   
										}
									}
									elsif($MinorGrooveSize > 3.264064)
									{
										print OUT "$sn\tPROMOTER SEQUENCE\n"; 
									}
								}
							}
							elsif($MinorGrooveSize > 3.290837)
							{
								if($MinorGrooveSize <= 3.295936)
								{
									if($AT <= 0.669323)
									{
										if($Flexibility_shift <= 5.820637)
										{
											if($Adeninecontent <= 0.577689)
											{
												print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
											}
											elsif($Adeninecontent > 0.577689)
											{
												if($Guaninecontent <= 0.318725)
												{
													print OUT "$sn\tPROMOTER SEQUENCE\n";  
												}
												elsif($Guaninecontent > 0.318725)
												{
													if($MinorGrooveSize <= 3.294821)
													{
														print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
													}
													elsif($MinorGrooveSize > 3.294821)
													{
														print OUT "$sn\tPROMOTER SEQUENCE\n"; 
													}
												}
											}
										}
										elsif($Flexibility_shift > 5.820637)
										{
											print OUT "$sn\tPROMOTER SEQUENCE\n";  
										}
									}
									elsif($AT > 0.669323)
									{
										print OUT "$sn\tPROMOTER SEQUENCE\n"; 
									}
								}
								elsif($MinorGrooveSize > 3.295936)
								{
									print OUT "$sn\tPROMOTER SEQUENCE\n";   
								}
							}
						}
					}
				}
				elsif($Rise_rise > 7.696737)
				{
					if($Six <= 293.383244)
					{
						if($acskew <= -0.043825)
						{
							if($Twist <= 36.63498)
							{
								print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
							}
							elsif($Twist > 36.63498)
							{
								if($acskew <= -0.10757)
								{
									print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
								}
								elsif($acskew > -0.10757)
								{
									if($Five <= 272.597886)
									{
										print OUT "$sn\tPROMOTER SEQUENCE\n";  
									}
									elsif($Five > 272.597886)
									{
										if($Twist <= 37.003625)
										{
											print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
										}
										elsif($Twist > 37.003625)
										{
											print OUT "$sn\tPROMOTER SEQUENCE\n";   
										}
									}
								}
							}
						}
						elsif($acskew > -0.043825)
						{
							if($FREEENERGY <= -303.18)
							{
								print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
							}
							elsif($FREEENERGY > -303.18)
							{
								if($Four <= 277.704042)
								{
									if($Mobilitymajorgroove <= 1.07498)
									{
										if($MinorGrooveSize <= 3.295936)
										{
											print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
										}
										elsif($MinorGrooveSize > 3.295936)
										{
											print OUT "$sn\tPROMOTER SEQUENCE\n";   
										}
									}
									elsif($Mobilitymajorgroove > 1.07498)
									{
										print OUT "$sn\tPROMOTER SEQUENCE\n";   
									}
								}
								elsif($Four > 277.704042)
								{
									print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
								}
							}
						}
					}
					elsif($Six > 293.383244)
					{
						if($Flexibility_shift <= 5.602948)
						{
							if($Mobilitymajorgroove <= 1.074263)
							{
								print OUT "$sn\tPROMOTER SEQUENCE\n";   
							}
							elsif($mobilitymajorgroove > 1.074263)
							{
								print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
							}
						}
						elsif($Flexibility_shift > 5.602948)
						{
							if($SCORE <= 4.542233)
							{
								print OUT "$sn\tPROMOTER SEQUENCE\n"; 
							}
							elsif($SCORE > 4.542233)
							{
								print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
							}
						}
					}
				}
			}
		}
		elsif($MinorGrooveSize > 3.29992)
		{
			if($AT <= 0.262948)
			{
				print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
			}
			elsif($AT > 0.262948)
			{
				if($acskew <= 0.139442)
				{
					if($Mobilitymajorgroove <= 1.050398)
					{
						if($Six <= 285.670171)
						{
							if($acskew <= 0.099602)
							{
								if($Direction <= 6.896414)
								{
									if($acskew <= 0.027888)
									{
										print OUT "$sn\tPROMOTER SEQUENCE\n"; 
									}
									elsif($acskew > 0.027888)
									{
										if($Four <= 262.570428)
										{
											print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
										}
										elsif($Four > 262.570428)
										{
											if($acskew <= 0.035857)
											{
												print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
											}
											elsif($acskew > 0.035857)
											{
												print OUT "$sn\tPROMOTER SEQUENCE\n";   
											}
										}
									}
								}
								elsif($Direction > 6.896414)
								{
									print OUT "$sn\tPROMOTER SEQUENCE\n";   
								}
							}
							elsif($acskew > 0.099602)
							{
								if($Adeninecontent <= 0.326693)
								{
									print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
								}
								elsif($Adeninecontent > 0.326693)
								{
									if($cgskew <= -0.044118)
									{
										print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
									}
									elsif($cgskew > -0.044118)
									{
										if($Two <= 239.112072)
										{
											if($Direction <= 12.661355)
											{
												print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
											}
											elsif($Direction > 12.661355)
											{
												print OUT "$sn\tPROMOTER SEQUENCE\n"; 
											}
										}
										elsif($Two > 239.112072)
										{
											print OUT "$sn\tPROMOTER SEQUENCE\n"; 
										}
									}
								}
							}
						}
						elsif($Six > 285.670171)
						{
							print OUT "$sn\tPROMOTER SEQUENCE\n";   
						}
					}
					elsif($Mobilitymajorgroove > 1.050398)
					{
						if($Flexibility_shift <= 6.234024)
						{
							if($FREEENERGY <= -315.1)
							{
								if($Six <= 283.72866)
								{
									if($acskew <= 0.043825)
									{
										if($Guaninecontent <= 0.334661)
										{
											if($Guaninecontent <= 0.278884)
											{
												print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
											}
											elsif($Guaninecontent > 0.278884)
											{
												if($Guaninecontent <= 0.302789)
												{
													print OUT "$sn\tPROMOTER SEQUENCE\n";   
												}
												elsif($Guaninecontent > 0.302789)
												{
													if($cgskew <= 0.130435)
													{
														print OUT "$sn\tPROMOTER SEQUENCE\n";   
													}
													elsif($cgskew > 0.130435)
													{
														if($AT <= 0.613546)
														{
															if($Twist <= 36.020518)
															{
																print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
															}
															elsif($Twist > 36.020518)
															{
																if($Six <= 282.472008)
																{
																	print OUT "$sn\tPROMOTER SEQUENCE\n";   
																}
																elsif($Six > 282.472008)
																{
																	print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
																}
															}
														}
														elsif($AT > 0.613546)
														{
															print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
														}
													}
												}
											}
										}
										elsif($Guaninecontent > 0.334661)
										{
											if($Twist <= 36.141275)
											{
												if($Adeninecontent <= 0.585657)
												{
													if($Flexibility_shift <= 6.126972)
													{
														print OUT "$sn\tPROMOTER SEQUENCE\n";   
													}
													elsif($Flexibility_shift > 6.12697)
													{
														print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
													}
												}
												elsif($Adeninecontent > 0.585657)
												{
													if($Rise_rise <= 7.593133)
													{
														if($Four <= 264.554522)
														{
															print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
														}
														elsif($Four > 264.554522)
														{
															print OUT "$sn\tPROMOTER SEQUENCE\n";   
														}
													}
													elsif($Rise_rise > 7.593133)
													{
														print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
													}
												}
											}
												elsif($Twist > 36.141275)
											{
												if($acskew <= -0.155378)
												{
													if($Rise_rise <= 7.650713)
													{
														print OUT "$sn\tPROMOTER SEQUENCE\n";   
													}
													elsif($Rise_rise > 7.650713)
													{
														print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
													}
												}
												elsif($acskew > -0.155378)
												{
													if($Five <= 272.137622)
													{
														print OUT "$sn\tPROMOTER SEQUENCE\n";    
													}
													elsif($Five > 272.137622)
													{
														if($Thyminecontent <= 0.697211)
														{
															if($Five <= 273.637532)
															{
																if($SCORE <= 5.336039)
																{
																	if($Six <= 279.144372)
																	{
																		print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
																	}
																	elsif($Six > 279.144372)
																	{
																		if($Direction <= 0.346614)
																		{
																			print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
																		}
																		elsif($Direction > 0.346614)
																		{
																			print OUT "$sn\tPROMOTER SEQUENCE\n";  
																		}
																	}
																}
																elsif($SCORE > 5.336039)
																{
																	print OUT "$sn\tPROMOTER SEQUENCE\n";    
																}
															}
															elsif($Five > 273.637532)
															{
																print OUT "$sn\tPROMOTER SEQUENCE\n";  
															}
														}
														elsif($Thyminecontent > 0.697211)
														{
															print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
														}
													}
												}
											}
										}
									}
									elsif($acskew > 0.043825)
									{
										if($Four <= 268.655941)
										{
											if($Flexibility_shift <= 6.087092)
											{
												print OUT "$sn\tPROMOTER SEQUENCE\n";  
											}
											elsif($Flexibility_shift > 6.087092)
											{
												if($Rise_rise <= 7.596414)
												{
													print OUT "$sn\tNOT A PROMOTER SEQUENCE\n"; 
												}
												elsif($Rise_rise > 7.596414)
												{
													print OUT "$sn\tPROMOTER SEQUENCE\n";    
												}
											}
										}
										elsif($Four > 268.655941)
										{
											if($Guaninecontent <= 0.414343)
											{
												print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
											}
											elsif($Guaninecontent > 0.414343)
											{
												if($Three <= 264.257308)
												{
													print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
												}
												elsif($Three > 264.257308)
												{
													if($Two <= 236.928608)
													{
														print OUT "$sn\tPROMOTER SEQUENCE\n";  
													}
													elsif($Two > 236.928608)
													{
															print OUT "$sn\tNOT A PROMOTER SEQUENCE\n"; 
													}
												}
											}
										}
									}
								}
								elsif($Six > 283.72866)
								{
									if($acskew <= 0.059761)	
									{
										if($Rise_rise <= 7.684018)
										{
											if($Adeninecontent <= 0.454183)
											{
												if($Mobilitymajorgroove <= 1.050797)
												{
													print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
												}
												elsif($Mobilitymajorgroove > 1.050797)
												{
													print OUT "$sn\tPROMOTER SEQUENCE\n";    
												}
											}
											elsif($Adeninecontent > 0.454183)
											{
												print OUT "$sn\tPROMOTER SEQUENCE\n";    
											}
										}
										elsif($Rise_rise > 7.684018)
										{
											if($Rise_rise <= 7.684385)
											{
												print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
											}
											elsif($Rise_rise > 7.684385)
											{
												if($FREEENERGY <= -316.71)
												{
													print OUT "$sn\tPROMOTER SEQUENCE\n";   
												}
												elsif($FREEENERGY > -316.71)
												{
													if($AT <= 0.629482)
													{
														if($Thyminecontent <= 0.62948)
														{
															print OUT "$sn\tPROMOTER SEQUENCE\n";     
														}
														elsif($Thyminecontent > 0.629482)
														{
															print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
														}
													}
													elsif($AT > 0.629482)
													{
														print OUT "$sn\tPROMOTER SEQUENCE\n";    
													}
												}
											}
										}
									}
									elsif($acskew > 0.059761)
									{
										if($Six <= 292.630939)
										{
											if($Four <= 275.667197)
											{
												if($Rise_rise <= 7.554105)
												{
													print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
												}
												elsif($Rise_rise > 7.554105)
												{
													print OUT "$sn\tPROMOTER SEQUENCE\n";    
												}
											}
											elsif($Four > 275.667197)
											{
												print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
											}
										}
										elsif($Six > 292.630939)
										{
											print OUT "$sn\tPROMOTER SEQUENCE\n";    
										}
									}
								}
							}
							elsif($FREEENERGY > -315.1)
							{
								if($AT <= 0.661355)
								{
									if($Six <= 291.581578)
									{
										if($Guaninecontent <= 0.414343)
										{
											if($Rise_rise <= 7.712955)
											{
												if($Two <= 199.506709)
												{
													if($acskew <= -0.163347)
													{
														if($Five <= 263.817523)
														{
															print OUT "$sn\tPROMOTER SEQUENCE\n";
															    
														}
														elsif($Five > 263.817523)
														{
															print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";
															  
														}
													}
													elsif($acskew > -0.163347)
													{
														print OUT "$sn\tPROMOTER SEQUENCE\n";
														    
													}
												}
												elsif($Two > 199.506709)
												{
													if($AT <= 0.61753)
													{
														if($FREEENERGY <= -312.85)
														{
															print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
														}
														elsif($FREEENERGY > -312.85)
														{
															print OUT "$sn\tPROMOTER SEQUENCE\n";   
														}
													}
													elsif($AT > 0.61753)
													{
														if($MinorGrooveSize <= 3.315857)
														{
															if($Flexibility_shift <= 5.751833)
															{
																if($AT <= 0.657371)
																{														
																	if($SCORE <= 3.650318)
																	{
																		print OUT "$sn\tPROMOTER SEQUENCE\n";   
																	}
																	elsif($SCORE > 3.650318)
																	{
																		if($Adeninecontent <= 0.577689)
																		{
																			print OUT "$sn\tPROMOTER SEQUENCE\n";    
																		}
																		elsif($Adeninecontent > 0.577689)
																		{
																			print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
																		}
																	}
																}
																elsif($AT > 0.657371)
																{
																	if($MinorGrooveSize <= 3.303904)
																	{
																		if($Thyminecontent <= 0.677291)
																		{
																			print OUT "$sn\tPROMOTER SEQUENCE\n";   
																		}
																		elsif($Thyminecontent > 0.677291)
																		{
																			print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
																		}
																	}
																	elsif($MinorGrooveSize > 3.303904)
																	{
																		print OUT "$sn\tPROMOTER SEQUENCE\n";   
																	}
																}
															}
															elsif($Flexibility_shift > 5.751833)
															{
																if($Guaninecontent <= 0.386454)
																{
																	if($cgskew <= -0.104762)
																	{
																		print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
																	}
																	elsif($cgskew > -0.104762)
																	{
																		print OUT "$sn\tPROMOTER SEQUENCE\n";  
																	}
																}
																elsif($Guaninecontent > 0.386454)
																{
																	print OUT "$sn\tPROMOTER SEQUENCE\n";   
																}
															}
														}
														elsif($MinorGrooveSize > 3.315857)
														{
															print OUT "$sn\tPROMOTER SEQUENCE\n";   
														}
													}
												}
											}
											elsif($Rise_rise > 7.712955)
											{
												if($AT <= 0.657371)
												{
													if($SCORE <= 4.617433)
													{
														if($Guaninecontent <= 0.342629)
														{
															print OUT "$sn\tPROMOTER SEQUENCE\n";    
														}
														elsif($Guaninecontent > 0.342629)
														{
															print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
														}
													}
													elsif($SCORE > 4.617433)
													{
														print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
													}
												}
												elsif($AT > 0.657371)
												{
													print OUT "$sn\tPROMOTER SEQUENCE\n";  
												}
											}
										}
										elsif($Guaninecontent > 0.414343)
										{
											if($Three <= 260.069669)
											{
												print OUT "$sn\tPROMOTER SEQUENCE\n";    
											}
											elsif($Three > 260.069669)
											{
                         						if($Flexibility_shift <= 5.846653)
												{
													if($acskew <= -0.043825)
													{
														if($AT <= 0.649402)
														{
															print OUT "$sn\tPROMOTER SEQUENCE\n";    
														}
														elsif($AT > 0.649402)
														{
															print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
														}
													}
													elsif($acskew > -0.043825)
													{
														print OUT "$sn\tPROMOTER SEQUENCE\n";  
													}
												}
												elsif($Flexibility_shift > 5.846653)
												{
													print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
												}
											}
										}
									}
									elsif($Six > 291.581578)
									{
										print OUT "$sn\tPROMOTER SEQUENCE\n";    
									}
								}
								elsif($AT > 0.661355)
								{
									print OUT "$sn\tPROMOTER SEQUENCE\n";   
								}
							}
						}
						elsif($Flexibility_shift > 6.234024)
						{
							if($acskew <= 0.035857)
							{
								if($AT <= 0.625498)
								{
									if($Rise_rise <= 7.579296)
									{
										if($Six <= 278.361719)
										{
											print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
										}
										elsif($Six > 278.361719)
										{
											print OUT "$sn\tPROMOTER SEQUENCE\n";  
										}
									}
									elsif($Rise_rise > 7.579296)
									{
										if($SCORE <= 2.725952)
										{
											if($acskew <= 0.011952)
											{
												print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
											}
											elsif($acskew > 0.01195)
											{
												print OUT "$sn\tPROMOTER SEQUENCE\n";    
											}
										}
										elsif($SCORE > 2.725952)
										{
											print OUT "$sn\tPROMOTER SEQUENCE\n";   
										}
									}
								}
								elsif($AT > 0.625498)
								{
									print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
								}
							}
							elsif($acskew > 0.035857)
							{
								if($Six <= 289.898989)
								{
									if($Thyminecontent <= 0.434263)
									{
										print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
									}
									elsif($Thyminecontent > 0.434263)
									{
										if($FREEENERGY <= -340.97)
										{
											if($Rise_rise <= 7.556656)
											{
												print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
											}
											elsif($Rise_rise > 7.556656)
											{
												print OUT "$sn\tPROMOTER SEQUENCE\n";  
											}
										}
										elsif($FREEENERGY > -340.97)
										{
											if($$Two <= 214.120963)
											{
												if($AT <= 0.581673)
												{
													print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
												}
												elsif($AT > 0.581673)
												{
													if($SCORE <= 7.374234)
													{
														print OUT "$sn\tPROMOTER SEQUENCE\n";   
													}
													elsif($SCORE > 7.374234)
													{
														print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
													}
												}
											}
											elsif($Two > 214.120963)
											{
												print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
											}
										}
									}
								}
								elsif($Six > 289.898989)
								{
									if($SCORE <= 4.858325)
									{
										print OUT "$sn\tPROMOTER SEQUENCE\n";  
									}
									elsif($SCORE > 4.858325)
									{
										if($Adeninecontent <= 0.669323)
										{
											if($Adeninecontent <= 0.585657)
											{
												if($Adeninecontent <= 0.569721)
												{
													print OUT "$sn\tPROMOTER SEQUENCE\n";  
												}
												elsif($Adeninecontent > 0.569721)
												{
													print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
												}
											}
											elsif($Adeninecontent > 0.585657)
											{
												print OUT "$sn\tPROMOTER SEQUENCE\n";  
											}
										}
										elsif($Adeninecontent > 0.669323)
										{
											print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
										}
									}
								}
							}
						}
					}
				}
				elsif($acskew > 0.139442)
				{
					if($Six <= 291.502144)
					{
						if($Direction <= 7.187251)
						{
							if($Rise_rise <= 7.60362)
							{
								print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
							}
							elsif($Rise_rise > 7.60362)
							{
								if($MinorGrooveSize <= 3.443347)
								{
									if($Rise_rise <= 7.608822)
									{
										if($SCORE <= 3.868625)
										{
											print OUT "$sn\tPROMOTER SEQUENCE\n";   
										}
										elsif($SCORE > 3.868625)
										{
											print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
										}
									}
									elsif($Rise_rise > 7.608822)
									{
										print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
									}
								}
								elsif($MinorGrooveSize > 3.443347)
								{
									print OUT "$sn\tPROMOTER SEQUENCE\n";   
								}
							}
						}
						elsif($Direction > 7.187251)
						{
							if($acskew <= 0.211155)
							{
								if($Adeninecontent <= 0.501992)
								{
									if($Adeninecontent <= 0.326693)
									{
										print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
									}
									elsif($Adeninecontent > 0.326693)
									{
										if($cgskew <= 0.316456)
										{
											print OUT "$sn\tPROMOTER SEQUENCE\n";    
										}
										elsif($cgskew > 0.316456)
										{
											print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";   
										}
									}
								}
								elsif($Adeninecontent > 0.501992)
								{
								if($Two <= 249.550153)
									{
										if($Flexibility_shift <= 6.059243)
										{
											print OUT "$sn\tPROMOTER SEQUENCE\n";   
										}
										elsif($Flexibility_shift > 6.059243)
										{
											print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
										}
									}
									elsif($Two > 249.550153)
									{
										print OUT "$sn\tPROMOTER SEQUENCE\n";    
									}
								}
							}
							elsif($acskew > 0.211155)
							{
								print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
							}
						}
					}
					elsif($Six > 291.502144)
					{
						if($acskew <= 0.227092)
						{
							if($Mobilitymajorgroove <= 1.05745)
							{
								if($cgskew <= 0.333333)
								{
									if($FREEENERGY <= -360.66)
									{
										print OUT "$sn\tPROMOTER SEQUENCE\n";  
									}
									elsif($FREEENERGY > -360.66)
									{
										if($Rise_rise <= 7.544586)
										{
											print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
										}
										elsif($Rise_rise > 7.544586)
										{
											if($Adeninecontent <= 0.541833)
											{
												if($Two <= 236.893989)
												{
													print OUT "$sn\tPROMOTER SEQUENCE\n";    
												}
												elsif($Two > 236.893989)	
												{
													print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
												}
											}
										   elsif($Adeninecontent > 0.541833)
											{
												print OUT "$sn\tPROMOTER SEQUENCE\n";  
											}
										}
									}
								}
								elsif($cgskew > 0.333333)
								{
									if($Four <= 281.551649)
									{
										print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
									}
									elsif($Four > 281.551649)
									{
										print OUT "$sn\tPROMOTER SEQUENCE\n";   
									}
								}
							}
							elsif($Mobilitymajorgroove > 1.05745)
							{
								if($Flexibility_shift <= 6.260677)
								{
									if($Direction <= -0.984064)
									{
										print OUT "$sn\tPROMOTER SEQUENCE\n";   
									}
									elsif($Direction > -0.984064)
									{
										if($Twist <= 35.980518)
										{
											print OUT "$sn\tPROMOTER SEQUENCE\n";    
										}
										elsif($Twist > 35.980518)
										{
											print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
										}
									}
								}
								elsif($Flexibility_shift > 6.260677)
								{
									print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
								}
							}
						}
						elsif($acskew > 0.227092)
						{
							print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
						}
					}
				}
			}
		}
	}
	elsif($Six > 302.744367)
	{
		if($AT <= 0.737052)
		{
			if($Twist <= 35.719841)
			{
				if($AT <= 0.258964)
				{
					print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
				}
				elsif($AT > 0.258964)
				{
					if($Six <= 309.153716)
					{
						if($Flexibility_shift <= 5.901833)
						{
							print OUT "$sn\tPROMOTER SEQUENCE\n";    
						}
						elsif($Flexibility_shift > 5.901833)
						{
							print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
						}
					}
					elsif($Six > 309.153716)
					{
						if($Rise_rise <= 7.634527)
						{
							if($Adeninecontent <= 0.298805)
							{
								if($Two <= 250.115461)
								{	
									print OUT "$sn\tPROMOTER SEQUENCE\n";  
								}
								elsif($Two > 250.115461)
								{
									print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
								}
							}
							elsif($Adeninecontent > 0.298805)
							{
								print OUT "$sn\tPROMOTER SEQUENCE\n";    
							}
						}
						elsif($Rise_rise > 7.634527)
						{
							print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
						}
					}
				}
			}
			elsif($Twist > 35.719841)
			{
				if($MinorGrooveSize <= 3.295936)
				{
					if($Six <= 331.988293)
					{
						if($MinorGrooveSize <= 3.26008)
						{
							if($MinorGrooveSize <= 3.244143)
							{
								if($SCORE <= 7.939091)
								{
									if($Two <= 208.464986)
									{
										if($Two <= 172.151819)
										{
											print OUT "$sn\tPROMOTER SEQUENCE\n";   
										}
										elsif($Two > 172.151819)
										{
											if($AT <= 0.721116)
											{
												print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
											}
											elsif($AT > 0.721116)
											{
												if($MinorGrooveSize <= 3.231076)
												{
													print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
												}
												elsif($MinorGrooveSize > 3.231076)
												{
													if($Four <= 292.454605)
													{
														if($Mobilitymajorgroove <= 1.076414)
														{
															if($Six <= 312.334098)
															{
																print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
															}
															elsif($Six > 312.334098)
															{
																print OUT "$sn\tPROMOTER SEQUENCE\n";  
															}
														}
														elsif($Mobilitymajorgroove > 1.076414)
														{
															print OUT "$sn\tPROMOTER SEQUENCE\n";  
														}
													}
													elsif($Four > 292.454605)
													{
														print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
													}
												}
											}
										}
									}
									elsif($Two > 208.464986)
									{
										if($MinorGrooveSize <= 3.228207)
										{
											print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
										}
										elsif($MinorGrooveSize > 3.228207)
										{
											if($Rise_rise <= 7.64436)
											{
												if($Six <= 317.752297)
												{
													print OUT "$sn\tNOT A PROMOTER SEQUENCE\n"; 
												}
												elsif($Six > 317.752297)
												{
													if($Direction <= -3.920319)
													{
														if($cgskew <= -0.049505)
														{
															print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
														}
														elsif($cgskew > -0.049505)
														{
															print OUT "$sn\tPROMOTER SEQUENCE\n";    
														}
													}
													elsif($Direction > -3.920319)
													{
														print OUT "$sn\tPROMOTER SEQUENCE\n";    
													}
												}
											}
											elsif($Rise_rise > 7.64436)
											{
												print OUT "$sn\tPROMOTER SEQUENCE\n";    
											}
										}
									}
								}
								elsif($SCORE > 7.939091)
								{
									if($Mobilitymajorgroove <= 1.078088)
									{
										if($FREEENERGY <= -286.17)
										{
											if($Adeninecontent <= 0.729084)
											{
												print OUT "$sn\tPROMOTER SEQUENCE\n";      
											}
											elsif($Adeninecontent > 0.729084)
											{
												print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
											}
										}
										elsif($FREEENERGY > -286.17)
										{
											print OUT "$sn\tPROMOTER SEQUENCE\n";  
										}
									}
									elsif($Mobilitymajorgroove > 1.078088)
									{
										print OUT "$sn\tPROMOTER SEQUENCE\n";  
									}
								}
							}
							elsif($MinorGrooveSize > 3.244143)
							{
								if($AT <= 0.705179)
								{
									if($SCORE <= 5.321008)
									{	
										print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
									}
									elsif($SCORE > 5.321008)
									{
										if($Six <= 308.765658)
										{
											print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
										}
										elsif($Six > 308.765658)
										{
											print OUT "$sn\tPROMOTER SEQUENCE\n";  
										}
									}
								}
								elsif($AT > 0.705179)
								{
									if($AT <= 0.717131)
									{
										if($MinorGrooveSize <= 3.256096)
										{
											if($AT <= 0.709163)
											{
												if($Mobilitymajorgroove <= 1.07741)
												{
													if($Two <= 217.749802)
													{
														print OUT "$sn\tPROMOTER SEQUENCE\n";    
													}
													elsif($Two > 217.749802)
													{
														print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
													}
												}
												elsif($Mobilitymajorgroove > 1.07741)
												{
													print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
												}
											}
											elsif($AT > 0.709163)
											{
												if($MinorGrooveSize <= 3.252112)
												{
													if($AT <= 0.713147)
													{
														if($FREEENERGY <= -293.38)
														{
															print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
														}
														elsif($FREEENERGY > -293.38)
														{
															print OUT "$sn\tPROMOTER SEQUENCE\n";   
														}
													}
													elsif($AT > 0.713147)
													{
														if($MinorGrooveSize <= 3.248127)
														{
															if($Mobilitymajorgroove <= 1.070956)
															{
																print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
															}
															elsif($Mobilitymajorgroove > 1.070956)
															{
																if($Three <= 281.395927)
																{
																
																	print OUT "$sn\tPROMOTER SEQUENCE\n";  
																}
																elsif($Three > 281.395927)
																{
																	print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
																}
															}
														}
														elsif($MinorGrooveSize > 3.248127)
														{
															print OUT "$sn\tPROMOTER SEQUENCE\n";  
														}
													}
												}
												elsif($MinorGrooveSize > 3.252112)
												{
													print OUT "$sn\tPROMOTER SEQUENCE\n";   
												}
											}
										}
										elsif($MinorGrooveSize > 3.256096)
										{
											print OUT "$sn\tPROMOTER SEQUENCE\n";   
										}
									}
									elsif($AT > 0.717131)
									{
										print OUT "$sn\tPROMOTER SEQUENCE\n";    
									}
								}
							}
						}
						elsif($MinorGrooveSize > 3.26008)
						{
							if($AT <= 0.701195)
							{
								if($AT <= 0.669323)
								{
									if($Rise_rise <= 7.646442)
									{
										print OUT "$sn\tPROMOTER SEQUENCE\n";   
									}
									elsif($Rise_rise > 7.646442)
									{
										if($Twist <= 36.56506)
										{
											if($SCORE <= 8.570741)
											{
												if($Two <= 202.4321)
												{
													print OUT "$sn\tPROMOTER SEQUENCE\n";  
												}
												elsif($Two > 202.4321)
												{
													print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
												}
											}
											elsif($SCORE > 8.570741)
											{
												print OUT "$sn\tPROMOTER SEQUENCE\n";  
											}
										}
										elsif($Twist > 36.56506)
										{
											print OUT "$sn\tPROMOTER SEQUENCE\n";    
										}
									}
								}
								elsif($AT > 0.669323)
								{
									if($MinorGrooveSize <= 3.291952)
									{
									if($AT <= 0.673307)
										{
											if($MinorGrooveSize <= 3.287968)
											{
												print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
											}
											elsif($MinorGrooveSize > 3.287968)
											{
												if($SCORE <= 5.854508)
												{
													if($Mobilitymajorgroove <= 1.07255)
													{
														if($Six <= 303.82956)
														{
															print OUT "$sn\tPROMOTER SEQUENCE\n";   
														}
														elsif($Six > 303.82956)
														{
															if($Rise_rise <= 7.699905)
															{
																print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
															}
															elsif($Rise_rise > 7.699905)
															{
																if($FREEENERGY <= -304.85)
																{
																	print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
																}
																elsif($FREEENERGY > -304.85)
																{
																	print OUT "$sn\tPROMOTER SEQUENCE\n";   
																}
															}
														}
													}
													elsif($Mobilitymajorgroove > 1.07255)
													{
														print OUT "$sn\tPROMOTER SEQUENCE\n";    
													}
												}
												elsif($SCORE > 5.854508)
												{
													if($cgskew <= -0.330097)
													{
														print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
													}
													elsif($cgskew > -0.330097)
													{
														print OUT "$sn\tPROMOTER SEQUENCE\n";    
													}
												}
											}
										}
										elsif($AT > 0.673307)
										{
											if($Six <= 311.883279)
											{
												if($MinorGrooveSize <= 3.272032)
												{
													if($AT <= 0.693227)
													{
														if($Direction <= 8.673307)
														{
															if($FREEENERGY <= -294.29)
															{
																print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
															}
															elsif($FREEENERGY > -294.29)
															{
																if($FREEENERGY <= -292.7)
																{
																	print OUT "$sn\tPROMOTER SEQUENCE\n";  
																}
																elsif($FREEENERGY > -292.7)
																{
																	print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
																}
															}
														}
														elsif($Direction > 8.673307)
														{
															print OUT "$sn\tPROMOTER SEQUENCE\n";  
														}
													}
													elsif($AT > 0.693227)
													{
														if($MinorGrooveSize <= 3.264064)
														{
															if($Mobilitymajorgroove <= 1.071514)
															{
																print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
															}
															elsif($Mobilitymajorgroove > 1.071514)
															{
																if($SCORE <= 8.690763)
																{
																	print OUT "$sn\tPROMOTER SEQUENCE\n";   
																}
																elsif($SCORE > 8.690763)
																{
																	print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
																}
															}
														}
														elsif($MinorGrooveSize > 3.264064)
														{
															if($AT <= 0.697211)
															{
																if($SCORE <= 5.899237)
																{
																	if($MinorGrooveSize <= 3.268048)
																	{
																		print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
																	}
																	elsif($MinorGrooveSize > 3.268048)
																	{
																		print OUT "$sn\tPROMOTER SEQUENCE\n";   
																	}
																}
																elsif($SCORE > 5.899237)
																{
																	print OUT "$sn\tPROMOTER SEQUENCE\n";    
																}
															}
															elsif($AT > 0.697211)
															{
																print OUT "$sn\tPROMOTER SEQUENCE\n";    
															}
														}
													}
												}
												elsif($MinorGrooveSize > 3.272032)
												{
													if($AT <= 0.685259)
													{
														if($MinorGrooveSize <= 3.287968)
														{
															if($AT <= 0.677291)
															{
																if($Rise_rise <= 7.669051)
																{
																	if($Twist <= 36.507888)
																	{
																		print OUT "$sn\tPROMOTER SEQUENCE\n";  
																	}
																	elsif($Twist > 36.507888)
																	{
																		print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
																	}
																}
																elsif($Rise_rise > 7.669051)
																{
																	print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
																}
															}
															elsif($AT > 0.677291)
															{
																if($MinorGrooveSize <= 3.28)
																{
																	if($Four <= 281.550067)
																	{
																		print OUT "$sn\tPROMOTER SEQUENCE\n";   
																	}
																	elsif($Four > 281.550067)
																	{
																		print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
																	}
																}
																elsif($MinorGrooveSize > 3.28)
																{
																	if($AT <= 0.681275)
																	{
																		if($MinorGrooveSize <= 3.283984)
																		{
																			if($Mobilitymajorgroove <= 1.070757)
																			{
																				print OUT "$sn\tNOT A PROMOTER SEQUENCE\n"; 
																			}
																			elsif($Mobilitymajorgroove > 1.070757)
																			{
																				print OUT "$sn\tPROMOTER SEQUENCE\n";  
																			}
																		}
																		elsif($MinorGrooveSize > 3.283984)
																		{
																			print OUT "$sn\tPROMOTER SEQUENCE\n";    
																		}
																	}
																	elsif($AT > 0.681275)
																	{
																		print OUT "$sn\tPROMOTER SEQUENCE\n";    
																	}
																}
															}
														}
														elsif($MinorGrooveSize > 3.287968)
														{
															print OUT "$sn\tPROMOTER SEQUENCE\n";  
														}
													}
													elsif($AT > 0.685259)
													{
														print OUT "$sn\tPROMOTER SEQUENCE\n";   
													}
												}
											}
											elsif($Six > 311.883279)
											{
												if($acskew <= 0.091837)
												{
													if($MinorGrooveSize <= 3.266932)
													{
														if($AT <= 0.697211)
														{
															if($Thyminecontent <= 0.693227)
															{
																if($Three <= 280.01919)
																{
																	print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
																}
																elsif($Three > 280.01919)
																{
																	print OUT "$sn\tPROMOTER SEQUENCE\n";    
																}
															}
															elsif($Thyminecontent > 0.693227)
															{
																print OUT "$sn\tPROMOTER SEQUENCE\n";    
															}
														}
														elsif($AT > 0.697211)
														{
															if($MinorGrooveSize <= 3.264064)
															{
																if($Rise_rise <= 7.683401)
																{
																	print OUT "$sn\tPROMOTER SEQUENCE\n";  
																}
																elsif($Rise_rise > 7.683401)
																{
																	print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
																}
															}
															elsif($MinorGrooveSize > 3.264064)
															{
																print OUT "$sn\tPROMOTER SEQUENCE\n";  
															}
														}
													}
													elsif($MinorGrooveSize > 3.266932)
													{
														if($MinorGrooveSize <= 3.286853)
														{
															if($AT <= 0.677291)
															{
																if($acskew <= 0.059761)
																{
																	print OUT "$sn\tPROMOTER SEQUENCE\n";  
																}
																elsif($acskew > 0.059761)
																{
																	print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
																}
															}
															elsif($AT > 0.677291)
															{
																if($Rise_rise <= 7.693472)
																{
																	print OUT "$sn\tPROMOTER SEQUENCE\n";  
																}
																elsif($Rise_rise > 7.693472)
																{
																	if($Two <= 201.923302)
																	{
																		print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
																	}
																	elsif($Two > 201.923302)
																	{
																		if($Rise_rise <= 7.693718)
																		{
																			print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
																		}
																		elsif($Rise_rise > 7.693718)
																		{
																			if($Five <= 313.77075)
																			{
																				print OUT "$sn\tPROMOTER SEQUENCE\n";    
																			}
																			elsif($Five > 313.77075)
																			{
																				print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
																			}
																		}
																	}
																}
															}
														}
														elsif($MinorGrooveSize > 3.286853)
														{
															print OUT "$sn\tPROMOTER SEQUENCE\n";  
														}
													}
												}
												elsif($acskew > 0.091837)
												{
													print OUT "$sn\tPROMOTER SEQUENCE\n";  
												}
											}
										}
									}
									elsif($MinorGrooveSize > 3.291952)
									{
										print OUT "$sn\tPROMOTER SEQUENCE\n";  
									}
								}
							}
							elsif($AT > 0.701195)
							{
								print OUT "$sn\tPROMOTER SEQUENCE\n";   
							}
						}
					}
					elsif($Six > 331.988293)
					{
						if($Five <= 346.769059)
						{
							if($MinorGrooveSize <= 3.248127)
							{
								if($AT <= 0.717131)
								{
									if($AT <= 0.713147)
									{
										print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
									}
									elsif($AT > 0.713147)
									{
										if($FREEENERGY <= -289.57)
										{
											if($Twist <= 36.525976)
											{
												print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
											}
											elsif($Twist > 36.525976)
											{
												if($Flexibility_shift <= 5.473904)
												{
													print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
												}
												elsif($Flexibility_shift > 5.473904)
												{
													print OUT "$sn\tPROMOTER SEQUENCE\n";    
												}
											}
										}
										elsif($FREEENERGY > -289.5)
										{
											print OUT "$sn\tPROMOTER SEQUENCE\n";    
										}
									}
								}
								elsif($AT > 0.717131)
								{
									if($MinorGrooveSize <= 3.232191)
									{
										if($AT <= 0.733068)
										{
											if($Adeninecontent <= 0.840637)
											{
												if($Guaninecontent <= 0.266932)
												{
													if($Flexibility_shift <= 5.239124)
													{
														print OUT "$sn\tPROMOTER SEQUENCE\n";    
													}
													elsif($Flexibility_shift > 5.239124)
													{
														print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
													}
												}
												elsif($Guaninecontent > 0.266932)
												{
													print OUT "$sn\tPROMOTER SEQUENCE\n";  
												}
											}
											elsif($Adeninecontent > 0.840637)
											{
												print OUT "$sn\tPROMOTER SEQUENCE\n";    
											}
										}
										elsif($AT > 0.733068)
										{
											print OUT "$sn\tPROMOTER SEQUENCE\n";    
										}
									}
									elsif($MinorGrooveSize > 3.232191)
									{
										if($MinorGrooveSize <= 3.244143)
										{
											if($AT <= 0.721116)
											{
												if($Six <= 341.948889)
												{
													if($Adeninecontent <= 0.641434)
													{
														print OUT "$sn\tPROMOTER SEQUENCE\n";    
													}
													elsif($Adeninecontent > 0.641434)
													{
														print OUT "$sn\tNOT A PROMOTER SEQUENCE\n"; 
													}
												}
												elsif($Six > 341.948889)
												{
													print OUT "$sn\tPROMOTER SEQUENCE\n";   
												}
											}
											elsif($AT > 0.721116)
											{
												print OUT "$sn\tPROMOTER SEQUENCE\n";   
											}
										}
										elsif($MinorGrooveSize > 3.244143)
										{
											print OUT "$sn\tPROMOTER SEQUENCE\n";   
										}
									}
								}
							}
							elsif($MinorGrooveSize > 3.248127)
							{
								print OUT "$sn\tPROMOTER SEQUENCE\n";  
							}
						}
						elsif($Five > 346.769059)
						{
							print OUT "$sn\tPROMOTER SEQUENCE\n";    
						}
					}
				}
				elsif($MinorGrooveSize > 3.295936)
				{
				
					if($acskew <= -0.112)
					{
						if($acskew <= -0.195219)
						{
							print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
						}
						elsif($acskew > -0.195219)
						{
							print OUT "$sn\tPROMOTER SEQUENCE\n";    
						}
					}
					elsif($acskew > -0.112)
					{
						if($Adeninecontent <= 0.661355)
						{
							if($acskew <= 0.322709)
							{
								print OUT "$sn\tPROMOTER SEQUENCE\n";   
							}
							elsif($acskew > 0.322709)
							{
								if($Guaninecontent <= 0.227092)
								{
									if($SCORE <= 4.324705)
									{
										print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
									}
									elsif($SCORE > 4.324705)
									{
										print OUT "$sn\tPROMOTER SEQUENCE\n";  
									}
								}
								elsif($Guaninecontent > 0.227092)
								{
									if($Three <= 284.763513)
									{
										if($Three <= 283.117868)
										{
											print OUT "$sn\tPROMOTER SEQUENCE\n";  
										}
										elsif($Three > 283.117868)
										{
											print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
										}
									}
									elsif($Three > 284.763513)
									{
										print OUT "$sn\tPROMOTER SEQUENCE\n";  
									}
								}
							}
						}
						elsif($Adeninecontent > 0.661355)
						{
							if($Rise_rise <= 7.534576)
							{
								if($Six <= 324.334511)
								{
									print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
								}
								elsif($Six > 324.334511)
								{
									if($Adeninecontent <= 0.752988)
									{
										print OUT "$sn\tPROMOTER SEQUENCE\n";   
									}
									elsif($Adeninecontent > 0.752988)
									{
										if($acskew <= 0.179283)
										{
											print OUT "$sn\tPROMOTER SEQUENCE\n";  
										}
										elsif($acskew > 0.179283)
										{
											print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
										}
									}
								}
							}
						   elsif($Rise_rise > 7.534576)
							{
								if($Six <= 327.196882)
								{
									if($Thyminecontent <= 0.414343)
									{
										if($acskew <= 0.211155)
										{
											if($Three <= 267.390989)
											{
												print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
											}
											elsif($Three > 267.390989)
											{
												print OUT "$sn\tPROMOTER SEQUENCE\n";   
											}
										}
										elsif($acskew > 0.211155)
										{
											print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
										}
									}
									elsif($Thyminecontent > 0.414343)
									{
										if($Flexibility_shift <= 5.541116)
										{
											if($Adeninecontent <= 0.808765)
											{
												if($Rise_rise <= 7.69836)
												{
													print OUT "$sn\tPROMOTER SEQUENCE\n";   
												}
												elsif($Rise_rise > 7.698368)
												{
													if($AT <= 0.657371)
													{
														if($Direction <= -2.115538)
														{
															print OUT "$sn\tPROMOTER SEQUENCE\n";  
														}
														elsif($Direction > -2.115538)
														{
															print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
														}
													}
													elsif($AT > 0.657371)
													{
														if($Adeninecontent <= 0.741036)
														{
															print OUT "$sn\tPROMOTER SEQUENCE\n";   
														}
														elsif($Adeninecontent > 0.741036)
														{	
															print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
														}
													}
												}
											}
											elsif($Adeninecontent > 0.808765)
											{
												print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
											}
										}
										elsif($Flexibility_shift > 5.541116)
										{
											if($acskew <= 0.155378)
											{
												print OUT "$sn\tPROMOTER SEQUENCE\n";  
											}
											elsif($acskew > 0.155378)
											{
												if($Rise_rise <= 7.647818)
												{
													if($cgskew <= 0.10679)
													{
														print OUT "$sn\tNOT A PROMOTER SEQUENCE\n"; 
													}
													elsif($cgskew > 0.106796)
													{
														if($Flexibility_shift <= 6.057371)
														{
															if($Rise_rise <= 7.626105)
															{
																if($Direction <= -5.788845)
																{
																	print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
																}
																elsif($Direction > -5.788845)
																{
																	if($SCORE <= 3.683055)
																	{
																		print OUT "$sn\tNOT A PROMOTER SEQUENCE\n"; 
																	}
																	elsif($SCORE > 3.683055)
																	{
																		print OUT "$sn\tPROMOTER SEQUENCE\n";     
																	}
																}
															}
															elsif($Rise_rise > 7.626105)
															{
																print OUT "$sn\tPROMOTER SEQUENCE\n";  
															}
														}
														elsif($Flexibility_shift > 6.057371)
														{
															print OUT "$sn\tPROMOTER SEQUENCE\n";  
														}
													}
												}
												elsif($Rise_rise > 7.647818)
												{
													print OUT "$sn\tPROMOTER SEQUENCE\n";    
												}
											}
										}
									}
								}
								elsif($Six > 327.196882)
								{
									print OUT "$sn\tPROMOTER SEQUENCE\n";   
								}
							}
						}
					}
				}
			}
		}
		elsif($AT > 0.737052)
		{
			if($Thyminecontent <= 0.717131)
			{
				if($Six <= 370.275412)
				{
					if($Direction <= -12.609562)
					{
						if($Three <= 295.4123)
						{
							print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
						}
						elsif($Three > 295.4123)
						{
							if($Mobilitymajorgroove <= 1.085139)
							{
								print OUT "$sn\tPROMOTER SEQUENCE\n";   
							}
							elsif($Mobilitymajorgroove > 1.085139)
							{
								if($Rise_rise <= 7.713015)
								{
									print OUT "$sn\tNOT A PROMOTER SEQUENCE\n"; 
								}
								elsif($Rise_rise > 7.713015)
								{
									print OUT "$sn\tPROMOTER SEQUENCE\n";  
								}
							}
						}
					}
					elsif($Direction > -12.609562)
					{
						if($acskew <= -0.016)
						{
							print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
						}
						elsif($acskew > -0.016)
						{
							print OUT "$sn\tPROMOTER SEQUENCE\n";  
						}
					}
				}
				elsif($Six > 370.275412)
				{
					if($MinorGrooveSize <= 3.199203)
					{
						if($Flexibility_shift <= 5.45239)
						{
							if($Thyminecontent <= 0.697211)
							{
								print OUT "$sn\tPROMOTER SEQUENCE\n";    
							}
							elsif($Thyminecontent > 0.697211)
							{
								if($Thyminecontent <= 0.701195)
								{
									print OUT "$sn\tNOT A PROMOTER SEQUENCE\n"; 
								}
								elsif($Thyminecontent > 0.701195)
								{
									print OUT "$sn\tPROMOTER SEQUENCE\n";    
								}
							}
						}
						elsif($Flexibility_shift > 5.45239)
						{
							if($SCORE <= 7.330424)
							{
								print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
							}
							elsif($SCORE > 7.330424)
							{
								print OUT "$sn\tPROMOTER SEQUENCE\n";  
							}
						}
					}
					elsif($MinorGrooveSize > 3.199203)
					{
						print OUT "$sn\tPROMOTER SEQUENCE\n";  
					}
				}
			}
			elsif($Thyminecontent > 0.717131)
			{
				if($MinorGrooveSize <= 3.212271)
				{
					if($AT <= 0.752988)
					{
						print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
					}
					elsif($AT > 0.752988)
					{
						if($MinorGrooveSize <= 3.207171)
						{
							if($cgskew <= 0.117647)
							{
								if($AT <= 0.760956)
								{
									print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
								}
								elsif($AT > 0.760956)
								{
									if($Mobilitymajorgroove <= 1.078725)
									{
										print OUT "$sn\tPROMOTER SEQUENCE\n";    
									}
									elsif($Mobilitymajorgroove > 1.078725)
									{
										print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
									}
								}
							}
							elsif($cgskew > 0.117647)
							{
								if($Rise_rise <= 7.59687)
								{
									print OUT "$sn\tPROMOTER SEQUENCE\n";   
								}
								elsif($Rise_rise > 7.59687)
								{
									print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
								}
							}
						}
						elsif($MinorGrooveSize > 3.207171)
						{
							if($Rise_rise <= 7.707026)
							{
								print OUT "$sn\tPROMOTER SEQUENCE\n";  
							}
							elsif($Rise_rise > 7.70702)
							{
								print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
							}
						}
					}
				}
				elsif($MinorGrooveSize > 3.212271)
				{
					if($MinorGrooveSize <= 3.224223)
					{
						if($AT <= 0.741036)
						{
							if($Direction <= -2.039841)
							{
								if($SCORE <= 9.38561)
								{
									print OUT "$sn\tPROMOTER SEQUENCE\n";   
								}
								elsif($SCORE > 9.38561)
								{
									print OUT "$sn\tNOT A PROMOTER SEQUENCE\n"; 
								}
							}
							elsif($Direction > -2.039841)
							{
								print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
							}
						}
						elsif($AT > 0.741036)
						{
							if($MinorGrooveSize <= 3.220239)
							{
								if($AT <= 0.74502)
								{
									if($FREEENERGY <= -285.8)
									{
										print OUT "$sn\tPROMOTER SEQUENCE\n";  
									}
									elsif($FREEENERGY > -285.8)
									{
										if($MinorGrooveSize <= 3.219124)
										{
											print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
										}
										elsif($MinorGrooveSize > 3.219124)
										{
											if($Adeninecontent <= 0.717131)
											{
												print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
											}
											elsif($Adeninecontent > 0.717131)
											{
												print OUT "$sn\tPROMOTER SEQUENCE\n";    
											}
										}
									}
								}
								elsif($AT > 0.74502)
								{
									if($AT <= 0.749004)
									{
										if($MinorGrooveSize <= 3.216255)
										{
											if($acskew <= 0.067729)
											{
												print OUT "$sn\tNOT A PROMOTER SEQUENCE\n";  
											}
											elsif($acskew > 0.067729)
											{
												print OUT "$sn\tPROMOTER SEQUENCE\n";    
											}
										}
										elsif($MinorGrooveSize > 3.216255)
										{
											print OUT "$sn\tPROMOTER SEQUENCE\n";   
										}
																			
									}
									elsif($AT > 0.74900)
									{
										print OUT "$sn\tPROMOTER SEQUENCE\n";  
									}
								}
								
							}							
							elsif($MinorGrooveSize > 3.220239)
							{
								print OUT "$sn\tPROMOTER SEQUENCE\n";     
							}
						}
					}
					elsif($MinorGrooveSize > 3.224223)
					{
						print OUT "$sn\tPROMOTER SEQUENCE\n";
						  
					}
				}
			}
		}
	}
}
}


close(IN);
close(OUT);