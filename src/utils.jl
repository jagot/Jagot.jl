function pdfmerge(output, args...)
    run(`pdfmerge $output $args`)
end

function ffmpeg_imgs(in_pattern, in_fps,
                     out_name, out_fps = 30,
                     vcodec = "h264", pix_fmt = "yuv420p",
                     assume_yes = true)
    run(`ffmpeg $(assume_yes ? "-y" : "") -framerate $in_fps -i $in_pattern -vcodec $vcodec -r $out_fps -pix_fmt $pix_fmt $out_name`)
end

export pdfmerge, ffmpeg_imgs
