function v = seeResults(field,filename, step)
  complete_name = sprintf('Results/%s/%s_%d.dat', filename, field, step)

  fid = fopen(complete_name, 'r');
  N = fread(fid,1,'int32');
  M = fread(fid,1,'int32');
  data = fread(fid,N*M,'double');
  v = zeros(M,N);
  for a = 1:1:N
    for b = 1:1:M
      v(b,a) = data((b-1) * N + a);
    end
  end
  figure1 = figure('visible','on');
  imagesc(v);
  set(gca,'YDir','normal');
  colorbar;
  fclose(fid);
end
