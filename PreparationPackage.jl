tpl = Template(;
    user="ChrisFuelOR",
    authors=["Christian Fuellner"],
    julia=v"1.5",
    plugins=[
        Git(; manifest=true),
        Codecov(),
        TravisCI(; x86=true),
        Documenter{TravisCI}(),
    ],
)

tpl("DynamicSDDiP")
