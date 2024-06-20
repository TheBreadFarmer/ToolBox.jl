function make_stackedhist(hists::Vector{Hist1D{T}};
    title::String="",
    xlabel::String="",
    ylabel::String="Entries (unweighted)",
    labels::Union{Nothing,Vector{String}}=nothing,
    xlimits=(nothing, nothing),
    ylimits=(0, nothing)
    ) where {T}

    set_theme!(theme_latexfonts())

    fig = Figure(; size=(1800,1200))

    ax = Axis(fig[1,1],
        title=title,
        titlesize=35,
        ylabel=ylabel,
        ylabelsize=20,
        xlabel=xlabel,
        xlabelsize=20
    )

    p = stackedhist!(hists; errors=false, error_color=:transparent)

    limits!(ax, xlimits[1], xlimits[2], ylimits[1], ylimits[2])

    ax.xticks = LinearTicks(10)
    ax.yticks = LinearTicks(7)

    if isnothing(labels)
        labels = ["Hist $(i)" for i in eachindex(hists)]
    else
        labels=labels
    end

    elements = [PolyElement(polycolor = p.attributes.color[][i]) for i in eachindex(labels)]

    Legend(fig[2,1], elements, labels; labelsize=20, orientation=:horizontal)

    return fig
end