Attribute VB_Name = "EDPEVO"
'
'  Calcul Numérique & Programmation
'  A/H.E. Lehtihet
'
'  Module   : << EDPEVO >>
'  Requis   :
'  Externes :
'
'  Descrip. : Equations aux dérivées partielles
'             Problèmes de conditions initiales
'
Option Explicit
Sub EDPEVO_1D(ByVal N As Integer, _
              ByRef U() As Double, _
              ByVal fac As Double, _
              ByVal M As Integer)

    ' Résolution de Ut = alpha*Uxx
    ' fac = alpha*dt/(dx*dx)
    
    Dim i As Integer, k As Integer
    Dim UC As Double, UG As Double, UD As Double
    Dim w As Double
    
    w = 1# - 2# * fac
    For k = 1 To M
        UG = U(0)
        UC = U(1)
        For i = 1 To N - 1
            UD = U(i + 1)
            U(i) = w * UC + fac * (UG + UD)
            UG = UC
            UC = UD
        Next i
    Next k
End Sub

