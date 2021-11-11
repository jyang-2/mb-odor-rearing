function data = read_raw_uint16(filepath, FOV)

    fid = fopen(filepath, 'r');
    if fid == -1
      error('Cannot open file: %s', FileName);
    end
    data = fread(fid, 'uint16');
    fclose(fid);
    data = uint16(data);
    data = reshape(data, FOV);
    data = pagetranspose(data);
end