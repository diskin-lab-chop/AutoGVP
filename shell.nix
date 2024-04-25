with import <nixpkgs> {};
{
  myProject = stdenv.mkDerivation {
    name = "AutoGVP";
    version = "1.0.0";
    src = if lib.inNixShell then null else nix;

    buildInputs = with rPackages; [
      bcftools
      R
      lubridate
      Biobase
      BiocManager
      optparse
      vcfR
      vroom
    ];
  };
}