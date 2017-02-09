using Glob

pgf_backend = matplotlib[:backends][:backend_pgf]

function add_icc(filename,
                 icc_prof = joinpath(dirname(@__FILE__), "sRGB_v4_ICC_preference.icc"))
    run(`convert $(filename) -profile $(icc_prof) $(filename)`)
end

function compile_pgf(filebase)
    w,h = [gcf()[s]() for s in
           [:get_figwidth,:get_figheight]]
    latex_source = """\\documentclass[12pt]{minimal}
\\usepackage[paperwidth=$(w)in, paperheight=$(h)in, margin=0in]{geometry}
\\usepackage{pdfx}
$(pgf_backend[:get_preamble]())
$(pgf_backend[:get_fontspec]())
\\usepackage{pgf}
\\begin{document}
\\centering
\\input{$(filebase).pgf}
\\end{document}"""
    open("$(filebase).tex", "w") do file
        write(file, latex_source)
    end
    tex_cmd = pgf_backend[:get_texcommand]()
    args = ["-interaction=nonstopmode", "-halt-on-error"]
    if tex_cmd in ["xetex", "xelatex"]
        println("WARNING: XeTeX is problematic")
        push!(args, "-shell-escape")
        push!(args, "-output-driver=xdvipdfmx -z 0")
    end
    run(`$(tex_cmd) $(args) $(filebase).tex`)
end

function save_pgf_with_icc(filename, args...;
                           icc_prof = joinpath(dirname(@__FILE__), "sRGB_v4_ICC_preference.icc"),
                           compile_pdf = true, kwargs...)
    savefig(filename, args...; kwargs...)
    dir,file,ext=dirname(filename),splitext(filename)...
    imgs = glob("$(file)-img*.png", dir)
    length(imgs) > 0 && println("Applying ICC profile to $(length(imgs)) image(s)")
    for img in imgs
        println(img)
        add_icc(img, icc_prof)
    end
    compile_pdf && compile_pgf(joinpath(dir, file))
end

export add_icc, save_pgf_with_icc
