function pdfmerge(output, args...)
    run(`pdfmerge $output $args`)
end

export pdfmerge
