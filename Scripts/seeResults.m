function v = seeResults(field,filename, step)
<<<<<<< HEAD
  complete_name = sprintf('/home/tom/Documents/Uliege/Master2/HPC/Project2/Results/%s/%s_%d.dat', filename, field, step)
=======
  complete_name = sprintf('Results/%s/%s_%d.dat', filename, field, step)
>>>>>>> master

  fid = fopen(complete_name, 'r');
  N = fread(fid,1,'int32');
  M = fread(fid,1,'int32');
<<<<<<< HEAD
  data = fread(fid, N*M,'double');
  v = zeros(M,N);
  for a = 1:1:N
    for b = 1:1:M
        v(b,a) = data((b-1) * N + a);
    end
  end
  figure('visible','on');
=======
  data = fread(fid,N*M,'double');
  v = zeros(M,N);
  for a = 1:1:N
    for b = 1:1:M
      v(b,a) = data((b-1) * N + a);
    end
    end1
  figure1 = figure('visible','on');
>>>>>>> master
  imagesc(v);
  set(gca,'YDir','normal');
  colorbar;
  fclose(fid);
end
