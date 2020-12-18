function[]=FUNC_SaveVideo(F,FR)
    n = 1;
    while true
        fn  ="Video"+string(n)+".avi";
        A = dir(char("**/"+fn));
        if size(A,1) ~= 0
            n = n+1;
        else
            break
        end
    end
    writerObj = VideoWriter(char(fn),'Motion JPEG AVI');
    writerObj.FrameRate = FR; 
    writerObj.Quality = 100; 
    open(writerObj);
    for i=1:length(F)
        frame = F(i) ;    
        writeVideo(writerObj, frame);
    end
    close(writerObj);
end