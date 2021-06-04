function calib_fileName= get_lower_calibFile(PicFile)

[DataDir, fileName]=fileparts(PicFile);
picNum= helper.getPicNum(fileName);

NH.calibFile= dir([DataDir filesep 'p*calib*raw*']);
if isempty(NH.calibFile)
    NH.calibFile= dir([DataDir filesep 'p*calib*']);
end

allCalibNums= cellfun(@(x) helper.getPicNum(x), {NH.calibFile.name}');

calibPicNum= find(allCalibNums<picNum, 1, 'last');
curDir= pwd;
cd(DataDir);
calib_fileName= [DataDir filesep helper.getFileName(allCalibNums(calibPicNum))];
cd(curDir);