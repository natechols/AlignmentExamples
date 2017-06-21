name := "nwalign"

scalaVersion in ThisBuild := "2.12.2"

version in ThisBuild := "0.1"

lazy val nwalign = (
  Project("nwalign", file("."))
    settings(
      libraryDependencies ++= Seq(
    //    "org.scala-lang.modules" %% "scala-xml" % "1.0.2",
        "org.scalatest" %% "scalatest" % "3.0.1"
      )
    ))
