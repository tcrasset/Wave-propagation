function makeMovie(field, filename, step, max)
  mov = VideoWriter('movie.avi');
  open(mov);
  for i = 2:step:max
    i
    complete_name = sprintf('Results/%s/%s__%d.dat', filename, field, i);
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
    figure1 = figure('visible','off');
    imagesc(v);
    set(gca,'YDir','normal');
    colorbar;
    frame = getframe(figure1);
    writeVideo(mov,frame);
  end
  close(mov);
end
