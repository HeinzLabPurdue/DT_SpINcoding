% File: startup.m
% Startup file for M. Heinz
% Created: 2005 Dec 23
%
% Modified for Purdue
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% List what projects to include in path
INCLUDE={'dummy'};  % Needed for book -keeping
INCLUDE{end+1}='Synology_HeinzLab';
INCLUDE{end+1}='JHU_Mfiles';
INCLUDE{end+1}='MHeinz_Laptop_Mfiles';

INCLUDE{end+1}='Recruitment';
INCLUDE{end+1}='LoudnessJNDsAnal';
INCLUDE{end+1}='RecruitPaper';




% Write out what was INCLUDED
disp(sprintf('\n**********\nINCLUDING:\n**********'))
disp(strvcat(INCLUDE{2:end}))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Setup path for each project included

if sum(strcmp('MHeinz_Laptop_Mfiles',INCLUDE))
    global Laptop_Mfiles_dir
    Laptop_Mfiles_dir='C:\HEINZ\Work\Research\MATLAB Mfiles';
    
    % Took out: June 03, 2007 - put whatever needed in MATLABfiles
    path(strcat(Laptop_Mfiles_dir,filesep,'EDYfiles'),path)
    path(strcat(Laptop_Mfiles_dir,filesep,'AFfunctions'),path)
    path(strcat(Laptop_Mfiles_dir,filesep,'Danilo - Data Import'),path)
    path(Laptop_Mfiles_dir,path)
end

if sum(strcmp('JHU_Mfiles',INCLUDE))
    global JHU_Mfiles_dir
    JHU_Mfiles_dir='C:\HEINZ\Work\JHU Cambridge PC\mgheinz\Recruitment\Mfiles';
    path(strcat(JHU_Mfiles_dir,filesep,'MGHfunctions'),path)
    path(JHU_Mfiles_dir,path)

end

if sum(strcmp('Synology_HeinzLab',INCLUDE))
    global HeinzLab_Mfiles_dir
    HeinzLab_Mfiles_dir='Y:\HeinzLabCode';
    path(strcat(HeinzLab_Mfiles_dir,filesep,'Behavior_Mfiles'),path)
    path(HeinzLab_Mfiles_dir,path)

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% early Purdue setup
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% INCLUDE{end+1}='NOHRstim';
%  INCLUDE{end+1}='R03anal';  % use 2.04
% INCLUDE{end+1}='STMPanal';  % use 2.04
% INCLUDE{end+1}='NManal';  
% INCLUDE{end+1}='NManal_new';  % Jun 17, 2008
% INCLUDE{end+1}='SLHS658';

if sum(strcmp('NOHRstim',INCLUDE))
   global NOHR_dir
   NOHR_dir='C:\Documents and Settings\Michael G. Heinz\My Documents\mgheinz\NOHR';
   eval(['cd ''' NOHR_dir filesep 'Stimuli'''])
   path(strcat(NOHR_dir,filesep,'Stimuli'),path)
   path(strcat(NOHR_dir,filesep,'Stimuli',filesep,'Vowel synthesis'),path)

   global FeaturesText FormsAtHarmonicsText InvertPolarityText
   FeaturesText={'F0','T0','F1','T1','F2','T2','F3','T3','TN'};
   FormsAtHarmonicsText={'yes','no'};
   InvertPolarityText={'no','yes'};
end

if sum(strcmp('R03anal',INCLUDE))
	%%% clean up NOHR/R03 names eventually, but left as is for now!!!!	
	global R03_dir
	R03_dir='C:\HEINZ\Work\Research\R03 Experiments';	
	global NOHR_dir NOHR_ExpList NOHR_IMPExpList
	eval(['cd ''' R03_dir filesep 'Data Analysis'''])
% 	path(strcat(R03_dir,filesep,'Data Analysis',filesep,'NOHR Mfiles'),path)
	path(strcat(R03_dir,filesep,'Data Analysis',filesep,'STMP Mfiles'),path)
	
	NOHR_dir=R03_dir;  % FIX THIS EVENTUALLY!!!
	NOHR_ExpList={'MH-2004_07_08-ANnorm-NOHR','MH-2004_08_02-AN-NOHR','MH-2004_09_02-ANnorm-NOHR','MH-2004_11_18-NOHRnorm', ...
		'MH-2004_12_14-NOHRdeafCat','MH-2005_02_10-ANmodel-ARO','MH-2005_03_28-ANnormal','MH-2005_04_18-ANnorm', ...
		'MH-2005_07_01-VCNnorm','MH-2005_07_08-VCNnorm','MH-2005_07_13-ANdeafcat'};
	NOHR_IMPExpList={'MH-2004_12_14-NOHRdeafCat','MH-2005_07_13-ANdeafcat'};

   global FeaturesText FormsAtHarmonicsText InvertPolarityText
   FeaturesText={'F0','T0','F1','T1','F2','T2','F3','T3','TN'};
   FormsAtHarmonicsText={'yes','no'};
   InvertPolarityText={'no','yes'};
	set(0,'DefaultTextInterpreter','tex')
end

if sum(strcmp('SLHS658',INCLUDE))
	% 	SLHS658_dir='C:\HEINZ\Work\Courses\SLHS 658 - Spring 2006\In-class Demos\Week1';
	SLHS658_dir='C:\HEINZ\Work\Courses\SLHS 658 - Spring 2006\SLHS 658 Toolbox';
	cd(SLHS658_dir)
	% 	ANmodel_dir='C:\HEINZ\Work\Courses\SLHS 658 - Spring 2006\In-class Demos\Week1\ARLO';
	% 	path(ANmodel_dir,path)
end

if sum(strcmp('STMPanal',INCLUDE))
	%%% clean up NOHR/R03 names eventually, but left as is for now!!!!	
	global STMP_dir
	STMP_dir='C:\HEINZ\Work\Research\R03 Experiments';	
	global STMP_ExpList STMP_IMPExpList 
	eval(['cd ''' STMP_dir filesep 'Data Analysis'''])
% 	path(strcat(STMP_dir,filesep,'Data Analysis',filesep,'NOHR Mfiles'),path)
	path(strcat(STMP_dir,filesep,'Data Analysis',filesep,'STMP Mfiles'),path)
	path(strcat(STMP_dir,filesep,'Data Analysis',filesep,'STMP Mfiles',filesep,'STMPfunctions'),path)
	path(strcat(STMP_dir,filesep,'Data Analysis',filesep,'STMP Mfiles',filesep,'SCCfunctions'),path)
	global MATLABfiles_dir
	path(strcat(MATLABfiles_dir,filesep,'NEL_TDT_Mfiles'),path)
	
	STMP_ExpList={'MH-2004_07_08-ANnorm-NOHR','MH-2004_08_02-AN-NOHR','MH-2004_09_02-ANnorm-NOHR','MH-2004_11_18-NOHRnorm', ...
		'MH-2004_12_14-NOHRdeafCat','MH-2005_02_10-ANmodel-ARO','MH-2005_03_28-ANnormal','MH-2005_04_18-ANnorm', ...
		'MH-2005_07_01-VCNnorm','MH-2005_07_08-VCNnorm','MH-2005_07_13-ANdeafcat','MH-2006_10_19-ANnorm', ...
      'MH-2006_11_03-ANnorm','MH-2006_11_16-ANnorm','MH-2007_02_23-ANnorm','MH-2007_03_09-ANnorm','MH-2007_04_13-ANnorm_wABR', ...
      'MH-2007_04_20-ANnorm_wABRs','AN-2007_06_20-modelSTMP','AN-2008_05_12-modelSTMP','AN-2008_05_19-modelBBN'};
   STMP_IMPExpList={'MH-2004_12_14-NOHRdeafCat','MH-2005_07_13-ANdeafcat'};

   global FeaturesText FormsAtHarmonicsText InvertPolarityText
   FeaturesText={'F0','T0','F1','T1','F2','T2','F3','T3','TN','NO'};
   FormsAtHarmonicsText={'yes','no'};
   InvertPolarityText={'no','yes'};
	set(0,'DefaultTextInterpreter','tex')
end

if sum(strcmp('NManal',INCLUDE))
	global NM_Exps NM_CALIBpics NM_Expdir NM_Analdir
	NM_Expdir='C:\HEINZ\Work\Research\R03 Experiments';	
	NM_Analdir='C:\HEINZ\Work\Papers\_PUBLISHED\2009\Neural Metrics - Ganesh\AN data';	
	
	NM_Exps={ ...
		'SK-2007_08_03-ANnorm', ...
		'SK-2007_09_13-AN_normal', ...
		'SK-2007_11_01-AN_normal', ...
		'SK-2007_08_30-ANexposed', ...
		'SK-2007_09_06-AN_impaired'};
	NM_CALIBpics=[2 9 9 3 6];

	path(NM_Analdir,path)
	eval(['cd ''' NM_Analdir ''''])
end

if sum(strcmp('NManal_new',INCLUDE))
	global ExpList NM_CALIBpics ROOT_dir InvertPolarityText
	ROOT_dir='C:\HEINZ\Work\Research\R03 Experiments';
	eval(['cd ''' ROOT_dir filesep 'Data Analysis'''])
	path(strcat(ROOT_dir,filesep,'Data Analysis',filesep,'STMP Mfiles'),path)
	path(strcat(ROOT_dir,filesep,'Data Analysis',filesep,'STMP Mfiles',filesep,'STMPfunctions'),path)
	path(strcat(ROOT_dir,filesep,'Data Analysis',filesep,'STMP Mfiles',filesep,'SCCfunctions'),path)
	global MATLABfiles_dir
	% 	path(strcat(MATLABfiles_dir,filesep,'NEL_TDT_Mfiles'),path)
	
	ExpList={ ...
		'SK-2007_08_03-ANnorm', ...
		'SK-2007_09_13-AN_normal', ...
		'SK-2007_11_01-AN_normal', ...
		'SK-2007_08_30-ANexposed', ...
		'SK-2007_09_06-AN_impaired', ...
		'MH-2008_06_03-ANnorm-chim', ...
		'SK-2008_12_04-AN_impaired'};
	NM_CALIBpics=[2 9 9 3 6 11];

% 	path('C:\HEINZ\Work\Papers\Neural Metrics - Ganesh\AN data',path)
% 	path('C:\HEINZ\Work\Research\R03 Experiments\Data Analysis\New Chin Analyses',path)
   InvertPolarityText={'no','yes'};
	set(0,'DefaultTextInterpreter','tex')
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% JHU setup
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% INCLUDE{end+1}='LoudnessJNDsAnal';
% INCLUDE{end+1}='RecruitPaper';

% INCLUDE{end+1}='CNexps';
% INCLUDE{end+1}='R03';
% INCLUDE{end+1}='JNDsPaper';
% INCLUDE{end+1}='ARO2004';
% INCLUDE{end+1}='IHCON2004';
% INCLUDE{end+1}='NoiseCorr';
% INCLUDE{end+1}='HY2004data';  % for Heinz and Young 2004 DATA

% INCLUDE{end+1}='Minnesota';
% INCLUDE{end+1}='Phase';
% INCLUDE{end+1}='RLFpaper';
% INCLUDE{end+1}='oldLoudness';
% INCLUDE{end+1}='ARO2003';
% INCLUDE{end+1}='ISH2003talk';
% INCLUDE{end+1}='ISH2003paper';

if sum(strcmp('HY2004data',INCLUDE))
	eval(['cd ''C:\HEINZ\Work\Research\Heinz and Young 2004 Data'''])
	
	disp('Not setting anything up for Heinz and Young 2004 data')
	
	%% Don't do anything!!!
end


%%%%%% AN Recruitment work
%%% DO FIRST, if any recruitment work is set to be done
if sum(strcmp('Recruitment',INCLUDE))
   global Recruit_dir 
	
%    Recruit_dir='C:\Documents and Settings\Michael G. Heinz\My Documents\mgheinz\Recruitment';
% 	Recruit_dir='C:\HEINZ\Work\Cambridge PC\mgheinz\Recruitment';
	Recruit_dir='C:\HEINZ\Work\JHU Cambridge PC\mgheinz\Recruitment';
	
	eval(['cd ''' Recruit_dir ''''])

   global SEVimpExps MILDimpExps NORMExps ALLimpExps
   SEVimpExps={'122001','032002','032702','040302'};
   MILDimpExps={'071002','091102','100902','101602','103002','110602'};
   NORMExps={'012501','013001','040501','041701','102501','110801','112901','061202'};
   ALLimpExps={'122001','032002','032702','040302','062002','071002','072302','091102','100902','101602', ...
         '103002','110602'};
   
   path(path,strcat(Recruit_dir,filesep,'POPData2'))
   path(path,strcat(Recruit_dir,filesep,'Mfiles'))
   path(path,strcat(Recruit_dir,filesep,'Mfiles',filesep,'MGHfunctions'))
   path(path,strcat(Recruit_dir,filesep,'Mfiles',filesep,'POPanal'))

   %   path(path,strcat(Recruit_dir,filesep,'Mfiles',filesep,'MGHfunctions',filesep,'Phase'))
end

if sum(strcmp('R03',INCLUDE))
   global R03_dir
   R03_dir='C:\Documents and Settings\Michael G. Heinz\My Documents\mgheinz\Grants\R03 - Apr 22, 2004\Mfiles';
   eval(['cd ''' R03_dir filesep ''''])
   path(R03_dir,path)
end

if sum(strcmp('RLFpaper',INCLUDE))
   global RLFpaper_dir
   RLFpaper_dir='C:\Documents and Settings\Michael G. Heinz\My Documents\mgheinz\Papers\RLFslopes';
   eval(['cd ''' RLFpaper_dir filesep 'Mfiles'''])
   path(strcat(RLFpaper_dir,filesep,'Mfiles'),path)
end

if sum(strcmp('LoudnessJNDsAnal',INCLUDE))
   global LoudnessJNDsAnal_dir
	%    LoudnessJNDsAnal_dir='C:\Documents and Settings\Michael G. Heinz\My Documents\mgheinz\Recruitment\Mfiles\LoudnessJND_MFiles';
	LoudnessJNDsAnal_dir='C:\HEINZ\Work\JHU Cambridge PC\mgheinz\Recruitment\Mfiles\LoudnessJND_MFiles';
	
	path(LoudnessJNDsAnal_dir,path)
   eval(['cd ''' LoudnessJNDsAnal_dir ''''])
end

if sum(strcmp('Phase',INCLUDE))
   path(strcat(Recruit_dir,filesep,'Mfiles',filesep,'Phase'),path)
   disp('IN Phase')
end

if sum(strcmp('JNDsPaper',INCLUDE))
   global JNDsPaper_dir
   JNDsPaper_dir='C:\Documents and Settings\Michael G. Heinz\My Documents\mgheinz\Papers\JNDs Paper';
   path(strcat(Recruit_dir,filesep,'Mfiles',filesep,'LoudnessJND_Mfiles'),path)
end

if sum(strcmp('RecruitPaper',INCLUDE))
   %% NEED TO INCLUDE: 'Recruitment','LoudnessJNDsAnal','RecruitPaper';
   if ~(sum(strcmp(INCLUDE,'Recruitment'))&sum(strcmp(INCLUDE,'LoudnessJNDsAnal'))&sum(strcmp(INCLUDE,'RecruitPaper')))
      error('Fix INCLUDES for RecruitPaper!!!!')
   end
     
   
   global RecruitPaper_dir
%    RecruitPaper_dir='C:\Documents and Settings\Michael G. Heinz\My Documents\mgheinz\Papers\RecruitmentAN Paper\Sub2';
   RecruitPaper_dir='C:\HEINZ\Work\JHU Cambridge PC\mgheinz\Papers\RecruitmentAN Paper - JARO 2005\Sub2';
   %    RecruitPaper_dir='C:\Documents and Settings\Michael G. Heinz\My Documents\mgheinz\Papers\RecruitmentAN Paper\Sub2\ExtraAnalysis\checkFlatHFloss';
   
   %    path(strcat(Recruit_dir,filesep,'Mfiles',filesep,'LoudnessJND_Mfiles'),path)
   path(strcat(RecruitPaper_dir,filesep,'Mfiles'),path)
   eval(['cd ''' RecruitPaper_dir filesep 'Mfiles'''])
end

if sum(strcmp('ARO2003',INCLUDE))
   global ARO2003_dir 
   ARO2003_dir='C:\Documents and Settings\Michael G. Heinz\My Documents\mgheinz\Presentations\ARO2003';
   eval(['cd ''' ARO2003_dir filesep 'Mfiles'''])
   path(strcat(ARO2003_dir,filesep,'Mfiles'),path)
end

if sum(strcmp('ARO2004',INCLUDE))
   global ARO2004_dir 
   ARO2004_dir='C:\Documents and Settings\Michael G. Heinz\My Documents\mgheinz\Presentations\ARO2004';
   eval(['cd ''' ARO2004_dir filesep 'Mfiles'''])
   path(strcat(ARO2004_dir,filesep,'Mfiles'),path)
end

if sum(strcmp('IHCON2004',INCLUDE))
   global IHCON2004_dir 
   IHCON2004_dir='C:\Documents and Settings\Michael G. Heinz\My Documents\mgheinz\Presentations\IHCON2004';
   eval(['cd ''' IHCON2004_dir filesep 'Talk' filesep 'Mfiles'''])
   path(strcat(IHCON2004_dir,filesep,'Talk',filesep,'Mfiles'),path)
   global LoudnessJNDsAnal_dir
   LoudnessJNDsAnal_dir='C:\Documents and Settings\Michael G. Heinz\My Documents\mgheinz\Recruitment\Mfiles\LoudnessJND_MFiles';
   path(LoudnessJNDsAnal_dir,path)
end

if sum(strcmp('ISH2003talk',INCLUDE))
   global ISHtalk_dir
   ISHtalk_dir='C:\Documents and Settings\Michael G. Heinz\My Documents\mgheinz\Papers\ISH2003\Talk';
   eval(['cd ''' ISHtalk_dir filesep 'Mfiles'''])
   path(strcat(ISHtalk_dir,filesep,'Mfiles'),path)
end

if sum(strcmp('ISH2003paper',INCLUDE))
   global ISHpaper_dir 
   ISHpaper_dir='C:\Documents and Settings\Michael G. Heinz\My Documents\mgheinz\Papers\ISH2003\Paper';
   path(strcat(ISHpaper_dir,filesep,'Mfiles'),path)
end

if sum(strcmp('oldLoudness',INCLUDE))
   oldLoudness_dir='C:\Documents and Settings\Michael G. Heinz\My Documents\mgheinz\Presentations\LabMtg_012403';
   path(strcat(oldLoudness_dir,filesep,'Mfiles'),path)
end

if sum(strcmp('NoiseCorr',INCLUDE))
   path(strcat(Recruit_dir,filesep,'Mfiles',filesep,'NoiseCorr'),path)
end

if sum(strcmp('Minnesota',INCLUDE))
   global Minn_dir
   Minn_dir='C:\Documents and Settings\Michael G. Heinz\My Documents\mgheinz\MinnesotaBME\Interview 051203\Seminar';
   path(strcat(Minn_dir,filesep,'Mfiles'),path)
end

if sum(strcmp('CNexps',INCLUDE))
   global CNexps_dir
   CNexps_dir='C:\Documents and Settings\Michael G. Heinz\My Documents\mgheinz\CNexps';
   eval(['cd ''' CNexps_dir filesep 'Data Analysis'''])
   path(strcat(CNexps_dir,filesep,'Data Analysis'),path)
   
   GE_dir='C:\Documents and Settings\Michael G. Heinz\My Documents\mgheinz\CNexps\NEL folders\Users\GE';
   GEMH_dir='C:\Documents and Settings\Michael G. Heinz\My Documents\mgheinz\CNexps\NEL folders\Users\GE_MH';
   path(GE_dir,path)
   path(GEMH_dir,path)
   
end
