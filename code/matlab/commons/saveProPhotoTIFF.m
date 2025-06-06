function saveProPhotoTIFF(rgbImg, outFile)
% saveProPhotoTiff  Save a ProPhoto-RGB image as 16-bit TIFF with embedded ICC
%
%   saveProPhotoTiff(rgbImg, outFile, iccPath)  
%     • rgbImg   : H×W×3 double in [0..1], already in ProPhoto RGB  
%     • outFile  : full path to write (e.g. 'plots/capture.tif')  

    % convert to 16-bit
    img16 = im2uint16(rgbImg);

    iccPath = '../../../data/icc/ProPhoto.icm';


    % read ICC profile bytes
    fid = fopen(iccPath,'r','ieee-be');
    if fid<0
        error("Cannot open ICC profile at %s", iccPath);
    end
    iccData = fread(fid, Inf, '*uint8');
    fclose(fid);

    % open Tiff and set tags
    t = Tiff(outFile,'w');
    tag.ImageLength         = size(img16,1);
    tag.ImageWidth          = size(img16,2);
    tag.Photometric         = Tiff.Photometric.RGB;
    tag.BitsPerSample       = 16;
    tag.SamplesPerPixel     = 3;
    tag.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
    tag.Compression         = Tiff.Compression.Deflate;         
    tag.RowsPerStrip        = min(512, size(img16,1));      % helps compression
    tag.ICCProfile          = iccData(:)';   % must be row uint8 array

    t.setTag(tag);
    t.write(img16);
    t.close();
end
